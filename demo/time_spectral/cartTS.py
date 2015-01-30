import os
import sys,getopt
import numpy as np
import time
sys.path.append('../../bin')
sys.path.append('../../../../../libTS/BUILD/libTS/bin')
from mpi4py import MPI
from cartInterface import cartModule
from libTSInterface import libTSModule
#
# Define Solver and Time-Spectral Modules
#
solver    = cartModule()
ts        = libTSModule()
#
# Initialize MPI
#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#
# Define Parameters
#
output_write  = 1
output_file   = 'stats.tmp'
nstep         = 500
ncheck        = 10
nsave         = nstep/ncheck
itern         = ncheck*np.arange(1,nsave+1) 
dq_res_two    = np.zeros(1)
rhs_res_two   = np.zeros(1)
dq_res_inf    = np.zeros(1)
rhs_res_inf   = np.zeros(1)
tcomp         = np.zeros(nsave)
qres          = np.zeros((nsave,2))
rres          = np.zeros((nsave,2))
#!
#> Preprocess
#!
#solver.inputParam(cfl,bc,jmax,kmax,lmax,irhs,ilhs)
solver.paramInput()
solver.initData()
ts.setData(solver.data)
#
#> Solve using AF scheme
#
s = 0
t = time.time()
for i in range(nstep):
    solver.rhs(viscous=False,bdf=False)
    ts.af_update()
    solver.lhs(i+1)
    if np.mod(i,ncheck) == 0:
        solver.computeNorm(rhs_res_two,dq_res_two,rhs_res_inf,dq_res_inf)
        tcomp[s] = time.time() - t
        rres[s,0]  = rhs_res_two
        qres[s,0]  = dq_res_two
        rres[s,1]  = rhs_res_inf
        qres[s,1]  = dq_res_inf
        s = s + 1
###############################################
#         Write Output File                   #
###############################################
t_ex = time.time()-t
t_pi = t_ex / nstep
if rank == 0:
    print "Execution time: ", t_ex," seconds"
    print "Iteration time: ", t_pi," seconds"
if rank == 0 and output_write == 1:
    print "Writing to file '", output_file,"'"
    file = open(output_file,"w")

    for i in range(s):
        file.write('%d\t%d\t%s\t%s\t%s\t%s\t%s\n' %(i+1,itern[i],tcomp[i],rres[i,0],rres[i,1],qres[i,0],qres[i,1]))
    file.close()
#
# Deallocate and finalize MPI
#
solver.outputData()
solver.finish()
ts.finish()

