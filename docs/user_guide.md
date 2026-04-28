# ase_fourstate User Guide

## Overview

`ase_fourstate` computes magnetic exchange coupling constants **J** using the **four-state method** with ASE (Atomic Simulation Environment) calculators. It supports any ASE-compatible calculator, and ships a custom `Abinit` calculator with automatic restart capabilities.

The exchange coupling is defined as:

```
J * S^2 = (E_FM_avg - E_AFM_avg) / 2
```

where four spin configurations (up-up, up-down, down-up, down-down) are evaluated on a pair of magnetic sites.

---

## Installation

```bash
pip install -e .
```

Requires Python >= 3.8, `ase`, and `numpy`.

---

## Quick Start

### One-command workflow with Abinit

```python
from ase_fourstate import run_abinit_exchange

results = run_abinit_exchange(
    structure_file="cri3_3x3.cif",
    magnetic_element="Cr",
    magnitude=3.0,
    max_neighbor=3,
    calc_params=dict(
        xc="PBE",
        ecut=1200,
        kpts=[4, 4, 1],
        pps="ONCV",
        smearing=["fermi-dirac", 0.001],
        toldfe=1e-6,
        nstep=100,
        diemac=4.0,
        nsppol=2,
    ),
    workdir="cri3_PBE",
)
```

This will:
1. Load the structure and set magnetic moments on `Cr` atoms.
2. Find neighbor shells up to `max_neighbor`.
3. Run Abinit calculations (N+2 instead of 4N).
4. Print and save a summary to `exchange_summary.txt`.

### Using any ASE calculator

```python
from ase.calculators.emt import EMT
from ase_fourstate import (
    load_structure,
    find_magnetic_site,
    collect_neighbors,
    compute_multi_neighbor_exchange,
    print_exchange_summary,
)

atoms = load_structure("structure.cif", "Cr", magnitude=3.0)
calc = EMT()
site1 = find_magnetic_site(atoms, "Cr")
neighbor_sites, neighbor_ids = collect_neighbors(site1, 3, "Cr", atoms)

results = compute_multi_neighbor_exchange(
    atoms, calc, site1, neighbor_sites, magnitude=3.0
)
print_exchange_summary(results, neighbor_ids)
```

---

## API Reference

### Structure and site helpers

| Function | Description |
|---|---|
| `load_structure(path, element, magnitude)` | Load CIF/POSCAR, set magnetic moments on `element` atoms |
| `find_magnetic_site(atoms, element)` | Return index of first atom of `element` |
| `find_neighbor(site1, neigh_idx, element, atoms)` | Find the `neigh_idx`-th neighbor shell atom (1-based) |
| `collect_neighbors(site1, max_n, element, atoms)` | Collect all neighbor shells up to `max_n`; returns `(sites, ids)` |
| `group_by_distance(dist_list, threshold)` | Group (index, distance) tuples into shells by proximity |

### Exchange calculation

#### `compute_four_state_exchange`

```python
compute_four_state_exchange(
    atoms, calculator, site1, site2,
    magnitude=1.0, verbose=True, neighbor_id="", use_nupdown=False
)
```

Performs the full 4-calculations (up-up, up-down, down-up, down-down) for a single pair. Returns a dict with `exchange_coupling` and the four energies.

#### `compute_multi_neighbor_exchange`

```python
compute_multi_neighbor_exchange(
    atoms, calculator, site1, neighbor_sites,
    magnitude=1.0, verbose=True, neighbor_ids=None, use_nupdown=False
)
```

Computes J for multiple neighbor shells in **N+2** calculations instead of 4N by exploiting symmetries:

- 1 reference (FM) calculation, shared by all pairs.
- 1 single-flip calculation, shared by all pairs.
- N double-flip calculations, one per neighbor.

Returns a dict keyed by neighbor label (e.g. `J1`, `J2`, ...), each containing `exchange_coupling`, `distance`, `site2`, and the three energies. A `_shared` key holds the reference and single-flip energies.

#### `print_exchange_summary`

```python
print_exchange_summary(results, neighbor_ids, structure_file="", outfile=None)
```

Prints a formatted table. If `outfile=True`, writes to `exchange_summary.txt`.

#### `run_abinit_exchange`

```python
run_abinit_exchange(
    structure_file, magnetic_element, magnitude, max_neighbor,
    calc_params, workdir=None, use_nupdown=False, dryrun=False
)
```

End-to-end workflow: load structure, create `Abinit` calculator, compute exchange, print/save summary. Handles directory management internally.

---

## The Abinit Calculator

The `ase_fourstate.Abinit` class is a custom `FileIOCalculator` for ABINIT with automatic restart support.

### Basic usage

```python
from ase_fourstate import Abinit

calc = Abinit(
    command="mpirun abinit < PREFIX.files > PREFIX.log",
    xc="PBE",
    ecut=1200,
    kpts=[4, 4, 1],
    pps="ONCV",
    toldfe=1e-6,
)
```

The command can also be set via the `ASE_ABINIT_SCRIPT` environment variable:

```bash
export ASE_ABINIT_SCRIPT="mpirun abinit < PREFIX.files > PREFIX.log"
```

Pseudopotential search paths are read from `ABINIT_PP_PATH` (colon-separated):

```bash
export ABINIT_PP_PATH="/path/to/pp/abinit"
```

### Automatic restart

When `reuse=True` (the default), the calculator:

1. Checks if a previous calculation converged. If so, reads results without re-running.
2. Checks for an existing wavefunction file (`o_WFK`). If found, symlinks it to `i_WFK` and sets `irdwfk=1`.
3. Falls back to the density file (`o_DEN`). If found, symlinks it to `i_DEN` and sets `irdden=1`.

Set `reuse=False` to always start from scratch.

### Supported pseudopotential families

| `pps` value | Description |
|---|---|
| `fhi` | FHI pseudopotentials |
| `hgh` | Hartwigsen-Goedecker-Hutter (lowest valence) |
| `hgh.sc` | HGH semicore |
| `hgh.k` | HGH with k-points |
| `tm` | Troullier-Martins |
| `paw` | PAW (hardest available) |
| `jth` | JTH PAW (standard) |
| `jth_sp` | JTH PAW (semicore/fine-core) |
| `gbrv` | GBRV |
| `ONCV` / `oncv` | ONCVPSP (optimal norm-conserving) |
| `GPAW` | GPAW PAW |
| `eric-gen` | Custom generated |

### XC functionals

Passed via the `xc` keyword. Maps to ABINIT `ixc` internally:

| `xc` | `ixc` |
|---|---|
| `LDA` | 7 |
| `PBE` | 11 |
| `revPBE` | 14 |
| `RPBE` | 15 |
| `WC` | 23 |
| `PBEsol` | -116133 |

Alternatively, pass `ixc` directly to bypass the mapping.

### Key parameters

| Parameter | Description |
|---|---|
| `xc` | Exchange-correlation functional |
| `ecut` | Plane-wave energy cutoff (eV) |
| `kpts` | Monkhorst-Pack grid, e.g. `[4, 4, 1]` |
| `gamma` | Gamma-point centered (default `True`) |
| `pps` | Pseudopotential family |
| `smearing` | Tuple `('fermi-dirac', width)` or `('gaussian', width)` |
| `toldfe` | Energy convergence tolerance (eV) |
| `nsppol` | Number of spin polarizations (2 for spin-polarized) |
| `nband` | Number of bands |
| `reuse` | Reuse converged results / restart from saved files (default `True`) |
| `charge` | System charge |
| `dryrun` | Write input files only, skip execution |

### DFT+U

```python
calc.set_Hubbard_U(
    {
        "Fe": {"L": 2, "U": 4.0, "J": 0.3},
        "O":  {"L": 1, "U": 1.0, "J": 0.3},
    },
    type=1,
)
```

- `L`: orbital angular momentum (0=s, 1=p, 2=d, 3=f)
- `U`: Hubbard U (eV)
- `J`: Hund's J (eV)
- `type`: LDA+U flavor (0=FLL without J, 1=FLL with J)

### Specialized calculations

```python
calc.relax_calculation(atoms, optcell=2, ntime=50)
calc.scf_calculation(atoms, dos=True)
calc.ldos_calculation(atoms, pdos=True)
calc.wannier_calculation(atoms, wannier_input)
```

---

## Tutorial: CrI3 with Abinit

```bash
python examples/run_cri3_abinit.py cri3_3x3.cif --cutoff 3 --ecut 1200 --dryrun
```

This creates the following directory structure, where each spin configuration gets its own subdirectory:

```
workdir/
  ref/          # FM reference state
  1flip/        # Single-flip state
  J1_2flip/     # Double-flip for J1
  J2_2flip/     # Double-flip for J2
  J3_2flip/     # Double-flip for J3
  psp -> /path/to/pseudopotentials
  exchange_summary.txt
```

Pseudopotential directories are automatically symlinked into each subdirectory.
