import os
import sys
sys.path.append('../../bin')
from cartInterface import cartModule
#
nsteps=10
nsubiter=5
#
solver=cartModule()
solver.paramInput()
solver.initData()
#
for i in range(nsteps):
    print "# time step :",i+1
    for p in range(1,nsubiter+1):
        solver.rhs()
        solver.lhs(p)
    solver.update_time()
solver.outputData()
solver.finish()
