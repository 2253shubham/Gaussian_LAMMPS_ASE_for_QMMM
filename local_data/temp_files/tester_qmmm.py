from ase import Atoms
from ase.io import read, write
from ase.calculators.qmmm import SimpleQMMM
from ase.calculators.gaussian import Gaussian
from ase.calculators.gromacs import Gromacs
import itertools as it

######################################################################################

def load_system(file, qm, mm):
    atoms = read(f'structures/{file}')
    CALC_QM = get_calculator(qm)
    CALC_MM = get_calculator(mm)
    return atoms, CALC_QM, CALC_MM
    # atoms.center(vacuum=4.0)


def get_calculator(calc_params, CALC_label):
    if 'gaussian' in CALC_label:
        # UPDATE SUB-TOPMON module add gaussian/g16.a03
        # prog='/share/apps/gaussian/g16.a03/g16'
        Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
        return Gaussian(label=CALC_label, **calc_params)
    elif 'lammmps' in CALC_label:
        LammpsCalculator.command = 'g16 < PREFIX.com > PREFIX.log'
        return LammpsCalculator(label=CALC_label, **calc_params)

   
def get_gaussian_params(**kwargs):
    params = dict(
        chk='gaussian.chk',
        method='b3lyp', basis='6-31G(d)', 
        # basisfile='gen.basis',
        # SCRF='SMD,Solvent=Water',
        # Geom='(Check, NewDefinition)',
        # Guess='Read',
        Integral='(Grid=UltraFine)', SCF='Tight,IntRep,XQC',
        **kwargs)
    # if 'nprocshared' not in params:
    #     params['nprocshared'] = sys.core_count
    # if 'mem' not in params and sys.mem_bytes is not None:
    #     params['mem'] = '{}GB'.format(int(sys.mem_bytes*0.6/(1024.**3)/(sys.core_count/params['nprocshared'])))
    return params


def get_lammps_params(**kwargs):
    params = dict()
    return params


######################################################################################
# define QM/MM calculators
QM_label = 'calc/gaussian'; MM_label = 'calc/lammps'
QM_params = get_gaussian_params(charge=0, mult=1)
MM_params = get_lammps_params(charge=0, mult=1)
# CALC_MM = Gromacs(doing_qmmm = True)
# CALC_MM.set_own_params_runs('extra_mdrun_parameters', ' -nt 1 ') # for serial runs

# Create system
file = 'mfi_T5/POSCAR'
atoms, CALC_QM, CALC_MM = load_system(file, [QM_params, QM_label], [MM_params, MM_label])
# olefin_wox = load_system('full_system.xsd')

# Make QM atoms selection with help of ase gui
qm_idx = [57, 61, 63, 88, 91, 110, 111, 112, 113]

######################################################################################
# Set up calculator: simple subtractive
atoms.SimpleQMMM(qm_idx, CALC_QM, CALC_MM)
atoms.get_potential_energy()
# olefin_wox.calc = SimpleQMMM(qm_idx, CALC_QM, CALC_MM)
# print(olefin_wox.get_potential_energy())
