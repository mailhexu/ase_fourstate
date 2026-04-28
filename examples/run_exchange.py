"""
Generic exchange profile calculation using any ASE calculator.

Usage:
    python run_exchange.py <structure_file> <element> <max_neighbor>

Example:
    python run_exchange.py cri3_3x3.cif Cr 3
"""

import sys
import os
from ase_fourstate import (
    load_structure,
    find_magnetic_site,
    collect_neighbors,
    compute_multi_neighbor_exchange,
    print_exchange_summary,
)


def compute_exchange_profile(
    structure_file,
    magnetic_element,
    max_neighbor,
    magnitude=3.0,
    calculator_name='emt',
    supercell=None,
):
    atoms = load_structure(structure_file, magnetic_element, magnitude)

    if supercell and supercell != [1, 1, 1]:
        atoms = atoms * supercell
        print(f"Created supercell {supercell}. Atoms: {len(atoms)}")

    calc = _make_calculator(calculator_name)
    if calc is None:
        return

    site1 = find_magnetic_site(atoms, magnetic_element)
    print(f"Reference magnetic site (site1): index {site1}, {atoms[site1].symbol}")

    print(f"\nIdentifying neighbor shells up to {max_neighbor}...")
    neighbor_sites, neighbor_ids = collect_neighbors(
        site1, max_neighbor, magnetic_element, atoms
    )

    if not neighbor_sites:
        print("No neighbor shells found. Exiting.")
        return

    N = len(neighbor_sites)
    print(f"\nComputing Exchange Parameters ({N + 2} calculations "
          f"instead of {4 * N})...")

    results = compute_multi_neighbor_exchange(
        atoms=atoms,
        calculator=calc,
        site1=site1,
        neighbor_sites=neighbor_sites,
        magnitude=magnitude,
        verbose=True,
        neighbor_ids=neighbor_ids,
    )

    print_exchange_summary(results, neighbor_ids, structure_file)


def _make_calculator(name):
    if name.lower() == 'emt':
        try:
            from ase.calculators.emt import EMT
            return EMT()
        except ImportError:
            print("EMT calculator not available.")
            return None
    print(f"Calculator '{name}' not implemented in this script.")
    return None


if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            struct = sys.argv[1]
            elem = sys.argv[2]
            cutoff = int(sys.argv[3])
        except IndexError:
            print("Usage: python run_exchange.py <structure_file> <element> <max_neighbor>")
            sys.exit(1)
    else:
        print("No arguments provided. Using default CrI3 settings for demonstration.")
        struct = 'cri3_3x3.cif'
        elem = 'Cr'
        cutoff = 3

    if os.path.exists(struct):
        compute_exchange_profile(
            structure_file=struct,
            magnetic_element=elem,
            max_neighbor=cutoff,
            calculator_name='emt',
        )
    else:
        print(f"Warning: {struct} not found.")
