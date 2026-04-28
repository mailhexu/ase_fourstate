
import os
import sys
from ase_fourstate import run_abinit_exchange

abinit_path = os.path.expanduser(
    "/home/ulg/phythema/hexu/.local/src/abinit-10.4.7/build/src/98_main"
)
os.environ['ASE_ABINIT_SCRIPT'] = f"mpirun {abinit_path}/abinit < PREFIX.files > PREFIX.log"
os.environ['ABINIT_PP_PATH'] = os.path.expanduser("~/.local/pp/abinit")


def run_mncl2_abinit_exchange(
    structure_file='MnCl2/mncl2_5x5_15A.cif',
    max_neighbor=3,
    magnitude=5.0,
    ecut=1400,
    kpts=[4, 4, 1],
    workdir='MnCl2',
    dryrun=False,
):
    """
    MnCl2 exchange with Abinit.

    Pseudopotentials (ONCVPSP-PW-PDv0.4):
        Mn: zion = 15  (3s2 3p6 3d5 4s2)
        Cl: zion =  7  (3s2 3p5)
    Valence electrons: 25*15 + 50*7 = 725
    nband = 464  (~28% above 725/2 = 363; 464 = 16*29, divisible by npband)
    """
    calc_params = dict(
        label='mncl2_abinit',
        xc='LDA',
        ecut=ecut,
        accuracy=3,
        optforces=0,
        kpts=kpts,
        pps='ONCV',
        smearing=['fermi-dirac', 0.001],
        toldfe=1e-6,
        nstep=100,
        diemac=4.0,
        nsppol=2,
        paral_kgb=1,
        np_spkpt=8,
        npfft=1,
        npband=16,
        nband=464,
        bandpp=1,
    )

    return run_abinit_exchange(
        structure_file=structure_file,
        magnetic_element='Mn',
        magnitude=magnitude,
        max_neighbor=max_neighbor,
        calc_params=calc_params,
        workdir=workdir,
        use_nupdown=True,
        dryrun=dryrun,
    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run Abinit Exchange Calculation for MnCl2")
    parser.add_argument('structure', nargs='?', default='MnCl2/mncl2_5x5_15A.cif',
                        help='Structure file')
    parser.add_argument('--cutoff', type=int, default=3, help='Max neighbor shell')
    parser.add_argument('--kpts', type=int, nargs=3, default=[4, 4, 1],
                        help='K-points (nkx nky nkz)')
    parser.add_argument('--ecut', type=float, default=1200, help='Energy cutoff (eV)')
    parser.add_argument('--dryrun', action='store_true', help='Write input files only, skip abinit execution')

    args = parser.parse_args()

    run_mncl2_abinit_exchange(
        structure_file=args.structure,
        max_neighbor=args.cutoff,
        kpts=args.kpts,
        ecut=args.ecut,
        dryrun=args.dryrun,
    )
