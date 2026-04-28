
import os
import sys
from ase_fourstate import run_abinit_exchange

abinit_path = os.path.expanduser(
    "/home/ulg/phythema/hexu/.local/src/abinit-10.4.7/build/src/98_main"
)
os.environ['ASE_ABINIT_SCRIPT'] = f"mpirun {abinit_path}/abinit < PREFIX.files > PREFIX.log"
os.environ['ABINIT_PP_PATH'] = os.path.expanduser("~/.local/pp/abinit")


def run_nii2_abinit_exchange(
    structure_file='NiI2/nii2_5x5_15A.cif',
    max_neighbor=3,
    magnitude=2.0,
    ecut=1400,
    kpts=[4, 4, 1],
    workdir='NiI2',
    dryrun=False,
):
    """
    NiI2 exchange with Abinit.

    Pseudopotentials (ONCVPSP-PW-PDv0.4):
        Ni: zion = 18  (3s2 3p6 3d8 4s2)
        I:  zion =  7  (5s2 5p5)
    Valence electrons: 25*18 + 50*7 = 800
    nband = 448  (~12% above 800/2 = 400; 448 = 16*28, divisible by npband)
    """
    calc_params = dict(
        label='nii2_abinit',
        ixc=7,
        ecut=ecut,
        accuracy=3,
        optforces=0,
        kpts=kpts,
        pps='ONCV',
        smearing=['fermi-dirac', 0.001],
        toldfe=1e-6,
        nstep=100,
        diemac=4.0,
        diemix=0.3,
        diemixmag=-4,
        nsppol=2,
        paral_kgb=1,
        np_spkpt=8,
        npfft=1,
        npband=16,
        nband=448,
        bandpp=1,
    )

    return run_abinit_exchange(
        structure_file=structure_file,
        magnetic_element='Ni',
        magnitude=magnitude,
        max_neighbor=max_neighbor,
        calc_params=calc_params,
        workdir=workdir,
        dryrun=dryrun,
    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run Abinit Exchange Calculation for NiI2")
    parser.add_argument('structure', nargs='?', default='NiI2/nii2_5x5_15A.cif',
                        help='Structure file')
    parser.add_argument('--cutoff', type=int, default=3, help='Max neighbor shell')
    parser.add_argument('--kpts', type=int, nargs=3, default=[4, 4, 1],
                        help='K-points (nkx nky nkz)')
    parser.add_argument('--ecut', type=float, default=1200, help='Energy cutoff (eV)')
    parser.add_argument('--dryrun', action='store_true', help='Write input files only, skip abinit execution')

    args = parser.parse_args()

    run_nii2_abinit_exchange(
        structure_file=args.structure,
        max_neighbor=args.cutoff,
        kpts=args.kpts,
        ecut=args.ecut,
        dryrun=args.dryrun,
    )
