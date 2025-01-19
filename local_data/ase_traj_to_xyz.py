from ase.io.trajectory import PickleTrajectory
import io, os
from ase.io import write

# reading in the trajectory file created during
# optimization
traj = PickleTrajectory("CO-Pd-ontop.traj")
# getting number of steps it took for geometry
# optimization
nsteps = dyn.get_number_of_steps()

string = "structure"

# get current working directory and make a scratch
# directory
path = os.getcwd()
path = path + "/scratch"
if not os.path.exists(path):
    os.makedirs(path)

# output file name
outFileName = "trajectory.xyz"
# write each structure from the .traj file in .xyz format
for i in range(0, nsteps + 1):
    atoms = traj[i]
    string = "structure%03d" % (i,) + ".xyz"
    outStruct = os.path.join(path, string)
    write(outStruct, atoms)
    # combines all optimization structures in one trajectory
    # file
    inFile = open(os.path.join(path, "structure%03d" % (i,) + ".xyz"), "r")
    fileStr = inFile.read()
    outFile = open(outFileName, "a")
    outFile.write(fileStr)
