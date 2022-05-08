# to import MPIManager
using MPIClusterManagers

# need to also import Distributed to use addprocs()
using Distributed

# specify, number of mpi workers, launch cmd, etc.
#manager=MPIWorkerManager(4)
manager=MPIWorkerManager(4)

# start mpi workers and add them as julia workers too.
addprocs(manager, exeflags="--project")