"""
ASE workflow to compute magnetic exchange coupling J using the 4-state method.
Reference: script/magnetic_exchange.py
"""

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.geometry import find_mic

def group_by_distance(dist_list, threshold=0.1):
    """
    Group a list of (index, distance) tuples by distance.
    
    Args:
        dist_list: list of tuples (idx, distance) or (idx, symbol, distance)
                   The function expects at least (idx, distance).
                   If the input tuples have 2 elements, it treats them as (idx, dist).
        threshold: maximum spread within one group (angstrom)

    Returns: 
        [
            {
                'average_dist': float,
                'items': [(idx, distance), ...]  # sorted by idx
            },
            ...
        ]
    """
    if not dist_list:
        return []

    # Handle different tuple lengths if necessary, but assuming (idx, dist) for simplicity
    # based on usage in find_neighbor.
    # The original script had (idx, symbol, distance) in docstring but (i, dist) in usage.
    # We will assume (idx, dist) is the minimal requirement and it's at the end.
    
    # Sort items by distance first (assuming distance is the last element)
    dist_list = sorted(dist_list, key=lambda x: x[-1])

    groups = []
    current_group = [dist_list[0]]

    for item in dist_list[1:]:
        if abs(item[-1] - current_group[0][-1]) <= threshold:
            # same group
            current_group.append(item)
        else:
            # close previous group
            groups.append(current_group)
            current_group = [item]

    # append final group
    groups.append(current_group)

    # Convert groups to list of dictionaries
    result = []
    for g in groups:
        avg = sum(x[-1] for x in g) / len(g)
        g_sorted = sorted(g, key=lambda x: x[0])  # sort by idx
        result.append({
            'average_dist': avg,
            'items': g_sorted
        })

    # Sort outer list by average distance
    result = sorted(result, key=lambda x: x['average_dist'])

    return result

def find_neighbor(site1, neigh_idx, filter_element, atoms):
    """
    Find the neighbor at a specific neighbor shell index.
    
    Args:
        site1: Index of the central atom.
        neigh_idx: Index of the neighbor shell (1-based). 1 means first nearest neighbor.
        filter_element: Symbol of the element to consider (e.g. 'Cr').
        atoms: ASE Atoms object.
        
    Returns:
        (site2, distance): Index of the neighbor atom and the distance.
    """
    positions = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()

    r0 = positions[site1]

    dist_list = []
    for i, r in enumerate(positions):
        if i == site1:
            continue

        if atoms[i].symbol != filter_element:
            continue

        # Minimum-image displacement vector
        dr = r - r0
        # find_mic returns a list of vectors for each input vector, here we have one input.
        # In some ASE versions it returns (vectors, length), in others just vectors.
        # find_mic signature: find_mic(D, cell, pbc) -> D_mic
        # D is list of vectors. D_mic is list of vectors.
        dr_mic = find_mic([dr], cell, pbc)[0][0]

        # Distance
        dist = np.linalg.norm(dr_mic)

        dist_list.append((i, dist))

    if not dist_list:
        raise ValueError(f"No neighbors found for element {filter_element}")

    # Sort by distance, then by index
    dist_list = sorted(dist_list, key=lambda x: (x[1], x[0]))
    dist_grouped_by_distance = group_by_distance(dist_list, threshold=0.1)

    if neigh_idx < 1 or neigh_idx > len(dist_grouped_by_distance):
        raise ValueError(f"Neighbor index {neigh_idx} out of range. Found {len(dist_grouped_by_distance)} shells.")

    group = dist_grouped_by_distance[neigh_idx-1]
    
    # Return the first atom in the sorted group
    return group['items'][0]

def compute_four_state_exchange(
    atoms: Atoms,
    calculator,
    site1: int,
    site2: int,
    magnitude: float = 1.0,
    verbose: bool = True,
    neighbor_id: str = "",
    use_nupdown: bool = False
) -> dict:
    """
    Compute the exchange coupling constant J*S^2 from four energy configurations.
    
    Uses the formula: J*S^2 = (E_FM_avg - E_AFM_avg) / 2
    where:
        E_FM_avg = (E_upup + E_downdown) / 2
        E_AFM_avg = (E_updown + E_downup) / 2
        
    This function performs 4 static calculations (no relaxation) corresponding to:
    up-up, up-down, down-up, down-down configurations on the two specified sites.
    
    Args:
        atoms: The ASE Atoms object (structure).
        calculator: The ASE calculator to use. Must be configured for spin-polarized calculations.
                    Note: If using file-based calculators (e.g. VASP, QE) in a script, 
                    ensure the calculator handles file IO correctly (e.g. overwriting files)
                    or manage directories externally if needed.
        site1: Index of the first magnetic site (0-based).
        site2: Index of the second magnetic site (0-based).
        magnitude: Magnitude of the magnetic moment (in Bohr magnetons).
        verbose: If True, print progress messages.
        use_nupdown: If True, set nupdown (total spin) on calculator for each configuration.
                     This constrains the total magnetic moment during the calculation.
        
    Returns:
        dict: A dictionary containing:
            - 'exchange_coupling': J*S^2 in eV
            - 'energy_upup': Energy of up-up configuration
            - 'energy_updown': Energy of up-down configuration
            - 'energy_downup': Energy of down-up configuration
            - 'energy_downdown': Energy of down-down configuration
    """
    
    # Validation
    num_sites = len(atoms)
    if site1 < 0 or site1 >= num_sites:
        raise ValueError(f"site1 index {site1} out of range [0, {num_sites-1}]")
    if site2 < 0 or site2 >= num_sites:
        raise ValueError(f"site2 index {site2} out of range [0, {num_sites-1}]")
        
    if verbose:
        print(f"Starting 4-state exchange calculation for sites {site1} and {site2}")
        print(f"Structure has {num_sites} atoms")
        print(f"Magnetization magnitude: {magnitude} µB")
    
    # Configurations
    # (sign_site1, sign_site2)
    configurations = {
        'upup': (1.0, 1.0),
        'updown': (1.0, -1.0),
        'downup': (-1.0, 1.0),
        'downdown': (-1.0, -1.0),
    }
    
    results = {}
    
    # Get initial magnetic moments from the input atoms to preserve other sites' moments
    # If not set, initialize to zeros
    try:
        initial_magmoms = atoms.get_initial_magnetic_moments()
    except (RuntimeError, AttributeError):
        initial_magmoms = np.zeros(num_sites)
    
    # Store original label/directory to restore later
    original_label = getattr(calculator, 'label', None)
    original_directory = getattr(calculator, 'directory', None)
    
    # If directory is '.', treat it as empty for cleaner subdirectory names
    if original_directory == '.':
        original_directory = None
    
    # Trigger psp directory creation by calling write_input once
    # This ensures psp exists before we create subdirectories
    import os
    import shutil
    
    # Save current directory
    saved_calc_directory = getattr(calculator, 'directory', None)
    
    # Temporarily use original directory to create psp
    if hasattr(calculator, 'directory'):
        calculator.directory = original_directory if original_directory else '.'
    
    # Write input once to trigger psp creation
    # We'll use a copy of atoms to avoid modifying the original
    try:
        temp_atoms = atoms.copy()
        temp_atoms.calc = calculator
        # Just write_input, don't calculate
        if hasattr(calculator, 'write_input'):
            calculator.write_input(temp_atoms)
    except Exception:
        pass  # Ignore errors, psp might already exist
    
    # Restore calculator directory
    if saved_calc_directory is not None and hasattr(calculator, 'directory'):
        calculator.directory = saved_calc_directory
    
    # Run calculations
    for config_name, (s1, s2) in configurations.items():
        # Prepare atoms for this configuration
        current_atoms = atoms.copy()
        
        # Set magnetic moments
        # We start with the base moments and modify only the target sites
        magmoms = initial_magmoms.copy()
        magmoms[site1] = s1 * magnitude
        magmoms[site2] = s2 * magnitude
        current_atoms.set_initial_magnetic_moments(magmoms)
        
        # Update the calculator's label for this configuration
        # This ensures each calculation runs in its own directory
        # Include neighbor_id if provided to distinguish J1, J2, etc.
        suffix = f"{neighbor_id}_{config_name}" if neighbor_id else config_name
        
        if hasattr(calculator, 'label') and original_label is not None:
            calculator.label = f"{original_label}_{suffix}"
        
        # Handle directory creation for each configuration
        if hasattr(calculator, 'directory'):
            import os
            import shutil
            
            # Determine the new subdirectory name
            if original_directory and original_directory != '.':
                new_directory = f"{original_directory}_{suffix}"
            else:
                # If no directory or '.', create subdirectory in current location
                new_directory = suffix
            
            calculator.directory = new_directory
            os.makedirs(new_directory, exist_ok=True)
            
            # Copy/symlink psp files
            # The psp files are typically in a 'psp' subdirectory of the current working dir
            cwd = os.getcwd()
            
            # Try multiple possible locations for psp directory
            possible_psp_paths = [
                os.path.join(cwd, 'psp'),  # psp in current directory
                'psp',  # relative path
            ]
            
            psp_src = None
            for path in possible_psp_paths:
                if os.path.exists(path):
                    psp_src = os.path.abspath(path)
                    break
            
            if psp_src:
                psp_dst = os.path.join(new_directory, 'psp')
                if not os.path.exists(psp_dst):
                    try:
                        # Create symlink (faster and saves space)
                        # Use relative path if possible for portability
                        rel_path = os.path.relpath(psp_src, os.path.abspath(new_directory))
                        os.symlink(rel_path, psp_dst)
                    except (OSError, NotImplementedError):
                        # Fall back to copying if symlink fails
                        shutil.copytree(psp_src, psp_dst)
        
        # Reset calculator to ensure it uses new directory
        if hasattr(calculator, 'reset'):
            calculator.reset()
        
        # Set nupdown if requested (fixed total magnetic moment)
        if use_nupdown:
            total_spin = np.sum(magmoms)
            if hasattr(calculator, 'set'):
                calculator.set(spinmagntarget=total_spin)
            elif hasattr(calculator, 'spinmagntarget'):
                calculator.spinmagntarget= total_spin
            if verbose:
                print(f"  Setting spinmagntarget = {total_spin:.1f}")
        
        # Attach calculator
        current_atoms.calc = calculator
        
        if verbose:
            print(f"Running configuration: {config_name} (site{site1}={magmoms[site1]:+.1f}, site{site2}={magmoms[site2]:+.1f})")
        
        # Double-check psp symlink exists before calculation
        # (The first iteration might not have it if psp dir wasn't created yet)
        if hasattr(calculator, 'directory'):
            import os
            import shutil
            cwd = os.getcwd()
            psp_src = os.path.join(cwd, 'psp')
            psp_dst = os.path.join(calculator.directory, 'psp')
            
            if os.path.exists(psp_src) and not os.path.exists(psp_dst):
                try:
                    rel_path = os.path.relpath(psp_src, os.path.abspath(calculator.directory))
                    os.symlink(rel_path, psp_dst)
                except (OSError, NotImplementedError):
                    shutil.copytree(psp_src, psp_dst)
            
        # Compute energy
        try:
            energy = current_atoms.get_potential_energy()
            results[f'energy_{config_name}'] = energy
            if verbose:
                print(f"  Energy: {energy:.6f} eV")
        except Exception as e:
            raise RuntimeError(f"Calculation failed for configuration {config_name}: {e}")
        
        # After each calculation, ensure psp symlinks exist for all directories
        # This fixes the case where the first calculation didn't have psp dir yet
        if hasattr(calculator, 'directory'):
            import os
            import shutil
            cwd = os.getcwd()
            psp_src = os.path.join(cwd, 'psp')
            
            # After first calculation, psp dir should exist in cwd
            # Create symlinks for any subdirs that are missing it
            if os.path.exists(psp_src):
                # Get all directories matching the neighbor pattern
                prefix = neighbor_id if neighbor_id else ''
                for item in os.listdir(cwd):
                    full_path = os.path.join(cwd, item)
                    if os.path.isdir(full_path) and (not prefix or item.startswith(prefix)):
                        psp_dst = os.path.join(full_path, 'psp')
                        if not os.path.exists(psp_dst):
                            try:
                                rel_path = os.path.relpath(psp_src, full_path)
                                os.symlink(rel_path, psp_dst)
                                if verbose:
                                    print(f"  Created psp symlink for {item}")
                            except (OSError, NotImplementedError):
                                shutil.copytree(psp_src, psp_dst)
                            except Exception:
                                pass  # Ignore errors, will retry next iteration
    
    # Restore original label/directory
    if original_label is not None:
        calculator.label = original_label
    # Restore directory (use '.' if it was None or '.')
    if hasattr(calculator, 'directory'):
        calculator.directory = original_directory if original_directory else '.'

    # Compute J*S^2
    e_upup = results['energy_upup']
    e_updown = results['energy_updown']
    e_downup = results['energy_downup']
    e_downdown = results['energy_downdown']
    
    e_fm_avg = (e_upup + e_downdown) / 2.0
    e_afm_avg = (e_updown + e_downup) / 2.0
    
    j_times_s2 = (e_fm_avg - e_afm_avg) / 2.0
    
    results['exchange_coupling'] = j_times_s2
    
    if verbose:
        print("-" * 40)
        print(f"Results:")
        print(f"  E_FM_avg : {e_fm_avg:.6f} eV")
        print(f"  E_AFM_avg: {e_afm_avg:.6f} eV")
        print(f"  J * S^2  : {j_times_s2:.6f} eV")
        print("-" * 40)
    
    return results

def _setup_calc_directory(calculator, suffix, original_label, original_directory, verbose=False):
    """
    Helper: set calculator label/directory for a given configuration suffix,
    create the subdirectory, and symlink psp files.
    
    Args:
        calculator: The ASE calculator.
        suffix: String suffix for the configuration (e.g. 'ref', '1flip', 'J1_2flip').
        original_label: Original calculator label to use as prefix.
        original_directory: Original calculator directory to use as prefix.
        verbose: Print progress messages.
    """
    import os
    import shutil

    if hasattr(calculator, 'label') and original_label is not None:
        calculator.label = f"{original_label}_{suffix}"

    if hasattr(calculator, 'directory'):
        if original_directory and original_directory != '.':
            new_directory = f"{original_directory}_{suffix}"
        else:
            new_directory = suffix

        calculator.directory = new_directory
        os.makedirs(new_directory, exist_ok=True)

        # Symlink psp files
        cwd = os.getcwd()
        psp_src = None
        for path in [os.path.join(cwd, 'psp'), 'psp']:
            if os.path.exists(path):
                psp_src = os.path.abspath(path)
                break

        if psp_src:
            psp_dst = os.path.join(new_directory, 'psp')
            if not os.path.exists(psp_dst):
                try:
                    rel_path = os.path.relpath(psp_src, os.path.abspath(new_directory))
                    os.symlink(rel_path, psp_dst)
                except (OSError, NotImplementedError):
                    shutil.copytree(psp_src, psp_dst)

    if hasattr(calculator, 'reset'):
        calculator.reset()


def _ensure_psp_symlinks(calculator, verbose=False):
    """
    Helper: ensure psp symlinks exist in the calculator's current directory.
    Called before running a calculation as a safety check.
    """
    import os
    import shutil

    if not hasattr(calculator, 'directory'):
        return

    cwd = os.getcwd()
    psp_src = os.path.join(cwd, 'psp')
    psp_dst = os.path.join(calculator.directory, 'psp')

    if os.path.exists(psp_src) and not os.path.exists(psp_dst):
        try:
            rel_path = os.path.relpath(psp_src, os.path.abspath(calculator.directory))
            os.symlink(rel_path, psp_dst)
        except (OSError, NotImplementedError):
            shutil.copytree(psp_src, psp_dst)


def _run_single_energy(atoms, calculator, magmoms, suffix,
                       original_label, original_directory, verbose=False,
                       use_nupdown=False):
    """
    Helper: run one static energy calculation with given magnetic moments.

    Args:
        atoms: Base ASE Atoms (will be copied).
        calculator: ASE calculator (label/directory will be modified in-place).
        magmoms: Array of magnetic moments for all atoms.
        suffix: Directory/label suffix for this configuration.
        original_label: Original calculator label.
        original_directory: Original calculator directory.
        verbose: Print progress.
        use_nupdown: If True, set nupdown on calculator based on total spin.

    Returns:
        float: Total energy in eV.
    """
    current_atoms = atoms.copy()
    current_atoms.set_initial_magnetic_moments(magmoms)

    _setup_calc_directory(calculator, suffix, original_label, original_directory, verbose)
    
    if use_nupdown:
        total_spin = np.sum(magmoms)
        if hasattr(calculator, 'set'):
            calculator.set(spinmagntarget=total_spin)
        elif hasattr(calculator, 'spinmagntarget'):
            calculator.spinmagntarget = total_spin
        if verbose:
            print(f"    Setting spinmagntarget = {total_spin:.1f}")
    
    current_atoms.calc = calculator
    _ensure_psp_symlinks(calculator, verbose)

    if verbose:
        print(f"  Running configuration: {suffix}")

    try:
        energy = current_atoms.get_potential_energy()
        if verbose:
            print(f"    Energy: {energy:.6f} eV")
        return energy
    except Exception as e:
        raise RuntimeError(f"Calculation failed for configuration {suffix}: {e}")


def _trigger_psp_creation(atoms, calculator, original_directory):
    """
    Helper: trigger psp directory creation by calling write_input once.
    """
    saved_calc_directory = getattr(calculator, 'directory', None)
    if hasattr(calculator, 'directory'):
        calculator.directory = original_directory if original_directory else '.'
    try:
        temp_atoms = atoms.copy()
        temp_atoms.calc = calculator
        if hasattr(calculator, 'write_input'):
            calculator.write_input(temp_atoms)
    except Exception:
        pass
    if saved_calc_directory is not None and hasattr(calculator, 'directory'):
        calculator.directory = saved_calc_directory


def compute_multi_neighbor_exchange(
    atoms: Atoms,
    calculator,
    site1: int,
    neighbor_sites: list,
    magnitude: float = 1.0,
    verbose: bool = True,
    neighbor_ids: list = None,
    use_nupdown: bool = False
) -> dict:
    """
    Compute exchange coupling J*S^2 for multiple neighbor shells efficiently.

    Exploits two symmetries of the four-state method to reduce 4N calculations
    to N+2:

    1. The reference (FM) state — all spins at +magnitude — is shared across
       all neighbor pairs.  Computed ONCE instead of N times.

    2. In a ferromagnetic reference background, flipping any single symmetry-
       equivalent magnetic atom gives the same total energy.  Therefore
       E(updown) = E(downup) for every pair, and these single-flip energies
       are all equal.  Computed ONCE instead of 2N times.

    Only the double-flip configuration (site1 AND site2 both flipped) is
    unique per neighbor pair, giving N calculations.

    Total: 1 (ref) + 1 (single-flip) + N (double-flips) = N + 2.

    The exchange coupling for neighbor n is:
        J_n * S^2 = (E_ref + E_2flip_n - 2 * E_1flip) / 4

    Args:
        atoms: ASE Atoms object with initial magnetic moments set for the
               FM reference state.
        calculator: ASE calculator (will be reused; label/directory are
                    modified between runs and restored afterwards).
        site1: Index of the central magnetic site (0-based).
        neighbor_sites: List of (site2_index, distance) tuples, one per
                        neighbor shell, as returned by find_neighbor().
        magnitude: Spin magnitude in µ_B.
        verbose: Print progress and energies.
        neighbor_ids: Optional list of string labels (e.g. ['J1','J2','J3']).
                      Defaults to ['J1', 'J2', ...].
        use_nupdown: If True, set nupdown (total spin) on calculator for each
                     configuration. This constrains the total magnetic moment.

    Returns:
        dict keyed by neighbor label, e.g.:
        {
            'J1': {
                'exchange_coupling': float,  # J*S^2 in eV
                'distance': float,
                'site2': int,
                'energy_ref': float,
                'energy_1flip': float,
                'energy_2flip': float,
            },
            ...
            '_shared': {
                'energy_ref': float,
                'energy_1flip': float,
                'num_calculations': int,  # = N+2
                'num_calculations_naive': int,  # = 4N (what would have been needed)
            }
        }
    """
    import os

    num_sites = len(atoms)
    N = len(neighbor_sites)

    if N == 0:
        raise ValueError("neighbor_sites must be non-empty")

    if neighbor_ids is None:
        neighbor_ids = [f"J{n+1}" for n in range(N)]
    if len(neighbor_ids) != N:
        raise ValueError("neighbor_ids length must match neighbor_sites length")

    # Validate sites
    if site1 < 0 or site1 >= num_sites:
        raise ValueError(f"site1 index {site1} out of range [0, {num_sites-1}]")
    for site2, dist in neighbor_sites:
        if site2 < 0 or site2 >= num_sites:
            raise ValueError(f"site2 index {site2} out of range [0, {num_sites-1}]")

    if verbose:
        print(f"Multi-neighbor exchange: site1={site1}, "
              f"{N} neighbors -> {N+2} calculations (naive: {4*N})")
        for nid, (s2, d) in zip(neighbor_ids, neighbor_sites):
            print(f"  {nid}: site2={s2}, dist={d:.3f} Å")

    # Get initial FM magnetic moments
    try:
        initial_magmoms = atoms.get_initial_magnetic_moments()
    except (RuntimeError, AttributeError):
        initial_magmoms = np.zeros(num_sites)

    original_label = getattr(calculator, 'label', None)
    original_directory = getattr(calculator, 'directory', None)
    if original_directory == '.':
        original_directory = None

    _trigger_psp_creation(atoms, calculator, original_directory)

    # ------------------------------------------------------------------
    # 1. Reference state: all spins at +magnitude (FM)
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n[1/{N+2}] Reference (FM) state")
    magmoms_ref = initial_magmoms.copy()
    magmoms_ref[site1] = magnitude  # ensure positive
    for site2, _ in neighbor_sites:
        magmoms_ref[site2] = magnitude
    energy_ref = _run_single_energy(
        atoms, calculator, magmoms_ref, "ref",
        original_label, original_directory, verbose, use_nupdown)

    # ------------------------------------------------------------------
    # 2. Single-flip state: flip site1 only (equivalent to flipping any
    #    single magnetic atom in the FM background)
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n[2/{N+2}] Single-flip state (site {site1} flipped)")
    magmoms_1flip = initial_magmoms.copy()
    magmoms_1flip[site1] = -magnitude
    for site2, _ in neighbor_sites:
        magmoms_1flip[site2] = magnitude
    energy_1flip = _run_single_energy(
        atoms, calculator, magmoms_1flip, "1flip",
        original_label, original_directory, verbose, use_nupdown)

    # ------------------------------------------------------------------
    # 3. Double-flip states: flip site1 AND site2(n) for each neighbor n
    # ------------------------------------------------------------------
    energy_2flip = {}
    for i, (nid, (site2, dist)) in enumerate(zip(neighbor_ids, neighbor_sites)):
        if verbose:
            print(f"\n[{3+i}/{N+2}] Double-flip state for {nid} "
                  f"(sites {site1} & {site2} flipped)")
        magmoms_2f = initial_magmoms.copy()
        magmoms_2f[site1] = -magnitude
        magmoms_2f[site2] = -magnitude
        energy_2flip[nid] = _run_single_energy(
            atoms, calculator, magmoms_2f, f"{nid}_2flip",
            original_label, original_directory, verbose, use_nupdown)

    # ------------------------------------------------------------------
    # Restore calculator state
    # ------------------------------------------------------------------
    if original_label is not None:
        calculator.label = original_label
    if hasattr(calculator, 'directory'):
        calculator.directory = original_directory if original_directory else '.'

    # ------------------------------------------------------------------
    # Extract J values
    # ------------------------------------------------------------------
    # J_n * S^2 = (E_ref + E_2flip_n - 2 * E_1flip) / 4
    #
    # Derivation:
    #   E_FM_avg  = (E_upup + E_downdown) / 2 = (E_ref + E_2flip_n) / 2
    #   E_AFM_avg = (E_updown + E_downup) / 2 = (E_1flip + E_1flip) / 2 = E_1flip
    #   J * S^2   = (E_FM_avg - E_AFM_avg) / 2
    #             = ((E_ref + E_2flip_n)/2 - E_1flip) / 2
    #             = (E_ref + E_2flip_n - 2*E_1flip) / 4
    # ------------------------------------------------------------------
    results = {}
    if verbose:
        print("\n" + "=" * 50)
        print("Results (multi-neighbor optimized)")
        print("=" * 50)
        print(f"  E_ref   = {energy_ref:.6f} eV")
        print(f"  E_1flip = {energy_1flip:.6f} eV")

    for nid, (site2, dist) in zip(neighbor_ids, neighbor_sites):
        e2 = energy_2flip[nid]
        j_s2 = (energy_ref + e2 - 2.0 * energy_1flip) / 4.0
        results[nid] = {
            'exchange_coupling': j_s2,
            'distance': dist,
            'site2': site2,
            'energy_ref': energy_ref,
            'energy_1flip': energy_1flip,
            'energy_2flip': e2,
        }
        if verbose:
            print(f"  {nid} (d={dist:.3f} Å): E_2flip={e2:.6f} eV  "
                  f"-> J*S^2 = {j_s2:.6f} eV")

    results['_shared'] = {
        'energy_ref': energy_ref,
        'energy_1flip': energy_1flip,
        'num_calculations': N + 2,
        'num_calculations_naive': 4 * N,
    }

    if verbose:
        savings = (1.0 - (N + 2) / (4 * N)) * 100
        print(f"\n  Calculations: {N+2} (saved {4*N - N - 2} "
              f"of {4*N}, {savings:.0f}% reduction)")
        print("=" * 50)

    return results


def load_structure(structure_file, magnetic_element, magnitude):
    """
    Load a structure file and set initial magnetic moments.

    Args:
        structure_file: Path to the structure file (CIF, POSCAR, etc.).
        magnetic_element: Symbol of the magnetic element (e.g. 'Cr').
        magnitude: Magnetic moment magnitude in Bohr magnetons.

    Returns:
        ASE Atoms object with initial magnetic moments set.
    """
    import os
    from ase.io import read as ase_read

    if not os.path.exists(structure_file):
        raise FileNotFoundError(f"Structure file '{structure_file}' not found.")

    print(f"Reading structure from {structure_file}...")
    atoms = ase_read(structure_file)
    elems = atoms.get_chemical_symbols()
    magmoms = [magnitude if elem == magnetic_element else 0.0 for elem in elems]
    atoms.set_initial_magnetic_moments(magmoms)
    print(f"Structure has {len(atoms)} atoms.")
    return atoms


def find_magnetic_site(atoms, magnetic_element):
    """
    Find the index of the first atom of the given magnetic element.

    Args:
        atoms: ASE Atoms object.
        magnetic_element: Symbol of the magnetic element.

    Returns:
        int: 0-based index of the first matching atom.

    Raises:
        ValueError: If no atom of the given element is found.
    """
    for i, atom in enumerate(atoms):
        if atom.symbol == magnetic_element:
            return i
    raise ValueError(f"No {magnetic_element} atoms found in structure.")


def collect_neighbors(site1, max_neighbor, magnetic_element, atoms):
    """
    Collect neighbor sites up to *max_neighbor* shells.

    Args:
        site1: Index of the central magnetic site.
        max_neighbor: Maximum neighbor shell to include (1-based).
        magnetic_element: Element symbol to filter neighbors.
        atoms: ASE Atoms object.

    Returns:
        (neighbor_sites, neighbor_ids)
        neighbor_sites: list of (site2_index, distance) tuples.
        neighbor_ids:   list of strings like ['J1', 'J2', ...].
    """
    neighbor_sites = []
    neighbor_ids = []
    for n in range(1, max_neighbor + 1):
        site2, dist = find_neighbor(site1, n, magnetic_element, atoms)
        neighbor_sites.append((site2, dist))
        neighbor_ids.append(f"J{n}")
        print(f"  Neighbor {n}: site={site2}, dist={dist:.3f} A")
    return neighbor_sites, neighbor_ids


def _format_exchange_summary(results, neighbor_ids, structure_file=""):
    """
    Build a formatted summary string of exchange coupling results.

    Args:
        results: Dict returned by compute_multi_neighbor_exchange.
        neighbor_ids: List of neighbor label strings (e.g. ['J1', 'J2']).
        structure_file: Optional structure file name for the header.

    Returns:
        str: The formatted summary text.
    """
    lines = []
    lines.append("=" * 60)
    if structure_file:
        lines.append(f"Results for {structure_file}")
    lines.append("=" * 60)

    shared = results['_shared']
    lines.append(f"E_ref   = {shared['energy_ref']:.6f} eV")
    lines.append(f"E_1flip = {shared['energy_1flip']:.6f} eV")
    lines.append("")

    lines.append(f"{'Shell':<6} {'Distance (A)':<14} {'E_2flip (eV)':<16} "
                 f"{'J*S^2 (eV)':<14}")
    lines.append("-" * 60)
    for nid in neighbor_ids:
        data = results[nid]
        lines.append(
            f"{nid:<6} {data['distance']:<14.4f} {data['energy_2flip']:<16.6f} "
            f"{data['exchange_coupling']:<14.6f}"
        )
    lines.append("-" * 60)
    lines.append(f"Total calculations: {shared['num_calculations']} "
                 f"(naive: {shared['num_calculations_naive']})")
    lines.append("=" * 60)
    return "\n".join(lines)


def print_exchange_summary(results, neighbor_ids, structure_file="",
                           outfile=None):
    """
    Print a formatted summary of exchange coupling results.

    Also writes the summary to *outfile* if given.

    Args:
        results: Dict returned by compute_multi_neighbor_exchange.
        neighbor_ids: List of neighbor label strings (e.g. ['J1', 'J2']).
        structure_file: Optional structure file name for the header.
        outfile: Optional path to write the summary to.
                 If ``True``, writes to ``exchange_summary.txt`` in the
                 current working directory.
    """
    import os

    text = _format_exchange_summary(results, neighbor_ids, structure_file)
    print("\n" + text)

    if outfile is True:
        outfile = "exchange_summary.txt"
    if outfile:
        with open(outfile, "w") as f:
            f.write(text + "\n")


def run_abinit_exchange(
    structure_file,
    magnetic_element,
    magnitude,
    max_neighbor,
    calc_params,
    workdir=None,
    use_nupdown=False,
    dryrun=False,
):
    """
    End-to-end workflow: load structure, create Abinit calculator, compute
    exchange couplings up to *max_neighbor* shells, and print a summary.

    Handles working-directory management (chdir / restore) internally.

    Args:
        structure_file: Path to the CIF/POSCAR file.
        magnetic_element: Symbol of the magnetic element.
        magnitude: Magnetic moment magnitude (µ_B).
        max_neighbor: Maximum neighbor shell index.
        calc_params: Dict of Abinit calculator parameters (xc, ecut, kpts, …).
                     Must NOT include ``command``; it is read from the
                     ``ASE_ABINIT_SCRIPT`` environment variable.
        workdir: Working directory for Abinit I/O.  Created if it does not
                 exist.  Defaults to the current directory.
        use_nupdown: Forward to compute_multi_neighbor_exchange.
        dryrun: If True, write input files but skip abinit execution.

    Returns:
        dict: Results from compute_multi_neighbor_exchange.
    """
    import os
    from .myabinit import Abinit

    atoms = load_structure(structure_file, magnetic_element, magnitude)

    cwd = os.getcwd()
    if workdir:
        abs_workdir = os.path.abspath(workdir)
        os.makedirs(abs_workdir, exist_ok=True)
        os.chdir(abs_workdir)

    try:
        cmd_str = os.environ['ASE_ABINIT_SCRIPT']
        calc = Abinit(command=cmd_str, dryrun=dryrun, **calc_params)
        print(f"Using command: {calc.command}")
        if dryrun:
            print("[DRYRUN] Abinit will NOT be executed; only input files will be written.")

        site1 = find_magnetic_site(atoms, magnetic_element)
        print(f"Computing J parameters for {magnetic_element} (Site {site1})...")

        neighbor_sites, neighbor_ids = collect_neighbors(
            site1, max_neighbor, magnetic_element, atoms
        )

        results = compute_multi_neighbor_exchange(
            atoms=atoms,
            calculator=calc,
            site1=site1,
            neighbor_sites=neighbor_sites,
            magnitude=magnitude,
            verbose=True,
            neighbor_ids=neighbor_ids,
            use_nupdown=use_nupdown,
        )

        print_exchange_summary(results, neighbor_ids, structure_file,
                               outfile=True)
        return results

    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    try:
        from ase.build import bulk
        from ase.calculators.emt import EMT

        print("Running test with EMT calculator on Fe BCC...")
        atoms = bulk('Fe', 'bcc', a=2.87) * (2, 1, 1)

        calc = EMT()

        res = compute_four_state_exchange(atoms, calc, site1=0, site2=1, magnitude=2.0)
        print("Test completed successfully.")

    except ImportError:
        print("ASE or EMT not available for testing.")
