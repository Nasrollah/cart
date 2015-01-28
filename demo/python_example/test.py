import os
import sys
import time
sys.path.append('../../bin')
from mpi4py import MPI
from cartInterface import cartModule
#
nsteps=10
nsubiter=40
#
solver=cartModule()
solver.paramInput()
solver.initData()
#
t = time.time()
for i in range(nsteps):
    if MPI.COMM_WORLD.Get_rank()==0:
	print "# time step :",i+1
    for p in range(1,nsubiter+1):
        solver.rhs()
        solver.lhs(p)
    solver.update_time()
if MPI.COMM_WORLD.Get_rank()==0:
    tcomp = time.time()-t
    print "Execution time: ",tcomp, " seconds"
solver.outputData()
solver.finish()

