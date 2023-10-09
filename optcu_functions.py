####################
from clusterx.structure import Structure

def total_energy(structure: Structure) -> float:
    """Compute total energy of optimized O @ PtCu system
    
    The structure is optimized using the BFGS method and the 
    total energy of the optimized structure is returned
    
    Parameters: 
    structure
        Structure for which the total energy is computed
        
    Returns:
        total energy
    """
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS
    from ase import Atoms
    from ase.io import read
    from ase.constraints import FixAtoms,FixedLine,FixedPlane
    from clusterx.utils import remove_vacancies
    import numpy as np

    nat = len(structure)
    at = structure.get_atoms()
    at.set_calculator(EMT())
    _structure = remove_vacancies(at)

    # Define constraints before relaxing:
    # - relax only top-most PtCu layer and O atoms
    # - restrict lateral displacements of O atoms
    # - restrict vertical displacements of top-most Pt and Cu atoms
    cstrs = []
    cstrs.append(FixAtoms(mask=[atom.z<11.0 for atom in _structure]))
    for i in np.argwhere(np.array([atom.z>15.0 for atom in _structure])==True):
        cstrs.append(FixedLine(i,direction=[0,0,1]))
    for i in np.argwhere(np.array([atom.z<15.0 and atom.z>13.0 for atom in _structure])==True):
        cstrs.append(FixedPlane(i,[0,0,1]))
    structure.set_constraint(cstrs)

    # Perform relaxation with BFGS and compute total energy
    e0 = _structure.get_potential_energy()
    rel = BFGS(_structure,trajectory = 'rel.traj',logfile='rel.log')
    rel.run(fmax=0.01)
    rel_str = read('rel.traj@-1')
    e1 = rel_str.get_potential_energy()

    return e1, rel_str, e0, structure
    #return e1
####################

####################
def references() -> None:
    """Compute reference energies
    
    Compute the reference energies needed for the 
    calculation of the energy of adsorption, namely,
    energy of O2 molecule (eO2), energy of Cu atom in bulk
    phase (eCu_bulk), energy of Pt atom in bulk phase (ePt_bulk),
    and energy of pristine Cu slab (eCu_slab).
    
    Return:
        eCu_bulk, ePt_bulk, eO2/2, eCu_slab, and string in reStructuredText (rst) format with information
    """
    from ase.build import bulk
    from ase.io.trajectory import Trajectory
    from ase.calculators.emt import EMT
    from ase import Atoms
    import numpy as np
    
    infostr = "\n==========================================\n"
    infostr +=  "Reference energies for the O @ PtCu system\n"
    infostr += "==========================================\n"
    
    infostr += "\n------\n"
    infostr += "  Pt  \n"
    infostr += "------\n\n"
    # Atomic Pt
    at = Atoms('Pt',[(0,0,0)])
    at.set_calculator(EMT())
    infostr += f"Energy Pt (atom): {at.get_potential_energy()}\n"

    # Bulk Pt
    at = bulk('Pt')
    at.set_calculator(EMT())
    cell = at.get_cell()
    traj = Trajectory('Pt.traj', 'w')
    for x in np.linspace(0.95, 1.05, 5):
        at.set_cell(cell * x, scale_atoms=True)
        at.get_potential_energy()
        traj.write(at)

    from ase.io import read
    from ase.eos import EquationOfState
    configs = read('Pt.traj@0:5')  # read 5 configurations
    # Extract volumes and energies:
    volumes = [at.get_volume() for at in configs]
    energies = [at.get_potential_energy() for at in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    infostr += f"Lattice constant Pt (bulk) = {np.power(v0*4,1/3)}\n"
    infostr += "\nEquation of state Pt bulk\n"
    infostr += "-------------------------\n\n"
    infostr += "+--------+--------+\n"
    infostr += "| Volume | Energy |\n"
    infostr += "+========+========+\n"
    for v, e in zip(volumes, energies):
        infostr += f"|{v:8.4f}|{e:8.4f}|\n"
    infostr += "+--------+--------+\n\n"
    
    infostr += f"Energy Pt (bulk): {e0}\n"
    ePt_bulk = e0

    infostr += "\n------\n"
    infostr += "  Cu  \n"
    infostr += "------\n\n"
    # Atomic Cu
    at = Atoms('Cu',[(0,0,0)])
    at.set_calculator(EMT())
    infostr += f"Energy Cu (atom): {at.get_potential_energy()}\n"

    # Bulk Cu
    at = bulk('Cu')
    at.set_calculator(EMT())
    cell = at.get_cell()
    traj = Trajectory('Cu.traj', 'w')
    for x in np.linspace(0.95, 1.05, 5):
        at.set_cell(cell * x, scale_atoms=True)
        at.get_potential_energy()
        traj.write(at)

    from ase.io import read
    from ase.eos import EquationOfState
    configs = read('Cu.traj@0:5')  # read 5 configurations
    # Extract volumes and energies:
    volumes = [at.get_volume() for at in configs]
    energies = [at.get_potential_energy() for at in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    infostr += f"Lattice constant Cu (bulk) = {np.power(v0*4,1/3)}\n"
    infostr += "\nEquation of state Cu bulk\n"
    infostr += "-------------------------\n\n"
    infostr += "+--------+--------+\n"
    infostr += "| Volume | Energy |\n"
    infostr += "+========+========+\n"
    for v, e in zip(volumes, energies):
        infostr += f"|{v:8.4f}|{e:8.4f}|\n"
    infostr += "+--------+--------+\n\n"
    
    infostr += f"Energy Cu (bulk): {e0}\n"
    eCu_bulk = e0

    infostr += "\n------\n"
    infostr += "  O   \n"
    infostr += "------\n\n"

    # Atomic O
    at = Atoms('O',[(0,0,0)])
    at.set_calculator(EMT())
    infostr += f"Energy O (atom): {at.get_potential_energy()}\n"

    # Molecular O2
    dmin=1.0
    dmax=1.2
    n=10
    es = []
    ds = []
    for i in range(n):
        d = dmin+i*(dmax-dmin)/(n-1)
        ds.append(d)
        o2 = Atoms('2O', [(0., 0., 0.), (0., 0., d)])
        o2.set_calculator(EMT())
        es.append(o2.get_potential_energy())
        
    infostr += "\nOptimization O2 molecule\n"
    infostr += "-------------------------\n\n"
    infostr += "+--------+--------+\n"
    infostr += "| Length | Energy |\n"
    infostr += "+========+========+\n"
    for d, e in zip(ds, es):
        infostr += f"|{d:8.4f}|{e:8.4f}|\n"
    infostr += "+--------+--------+\n\n"

    eO2 = min(es)
    infostr += f"Energy O2: {eO2}\n"
    infostr += f"Energy O: {eO2/2.0}"
    eCu_slab = 0.0
    return eCu_bulk, ePt_bulk, eO2/2.0, eCu_slab, infostr
####################

####################
from clusterx.structure import Structure

def ads_energy(structure, eCu_bulk=-6.978e-3, ePt_bulk=-1.821e-4, eO=3.134e-1, eCu_slab=11.2951):
    """Compute energy of adsorption of O @ PtCu surface alloy
    
    Parameters:
    structure
        The energy of adsorption of structure is computed
        
    Returns:
        Energy of adsorption per surface site
    """
    e,_,_,_ = total_energy(structure)
    nsites = structure.get_index()
    conc = structure.get_fractional_concentrations()

    xO = conc[0][1] # sublattice type 0: ['X' 'O'], concentration of Oxygen 
    xPt = conc[2][1] # sublattice type 2: ['Cu' 'Pt'], concentration of Pt

    eads = e/nsites - xO*eO - xPt*(ePt_bulk-eCu_bulk) - eCu_slab/nsites

    return eads
####################
