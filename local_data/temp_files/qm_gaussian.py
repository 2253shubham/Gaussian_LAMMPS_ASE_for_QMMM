from ase import Atoms
from ase.io import read, write
from ase.calculators.gaussian import Gaussian

######################################################################################

def load_system(file, calc_params, calc_label):
    #atoms = read(f'structures/{file}')
    atoms = read(f'structures/{file}.pdb')
    atoms.set_calculator(get_calculator(calc_params, calc_label))
    # atoms.center(vacuum=4.0)
    return atoms

def get_calculator(calc_params, calc_label):
    Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
    return Gaussian(label=calc_label, **calc_params)

def get_gaussian_params(**kwargs):
    params = dict(
    	mem='1000MW',
        #chk='gaussian.chk',
        #method='b3lyp', basis='6-31G(d)', 
        basisfile='gen_temp.basis',
        # SCRF='SMD,Solvent=Water',
        # Geom='(Check, NewDefinition)',
        # Guess='Read',
        #Integral='(Grid=UltraFine)', SCF='Tight,IntRep,XQC',
        #label='calc/gaussian',
        #xc='wB97X',
        #basis='6-31+G*',
        #scf='maxcycle=200',
        xc="M062X",
        basis="Gen/Auto",
        SCRF="SMD, Solvent=n-Octanol",
        pop='Hirshfeld',
        #Symmetry="None",
        **kwargs)
    # if 'nprocshared' not in params:
    #     params['nprocshared'] = sys.core_count
    # if 'mem' not in params and sys.mem_bytes is not None:
    #     params['mem'] = '{}GB'.format(int(sys.mem_bytes*0.6/(1024.**3)/(sys.core_count/params['nprocshared'])))
    return params
######################################################################################
# Create system

#file = 'mfi_T5/POSCAR'
file = 'MFIt23_with_ethane_QM_part'
calc_label = 'calc_temp/gaussian_stdaloneQMpart_temp'
calc_params = get_gaussian_params(charge=0, mult=1)
atoms = load_system(file, calc_params, calc_label)
atoms.get_potential_energy()
# print(atoms)