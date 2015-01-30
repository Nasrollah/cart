#
#
#               cart Object interface for HI-ARMS
#               created by J. Lefell and J. Sitaraman, 
#               last updated 12/1/2014
#
# Standard Python modules

import sys
import os
import copy
import string
import types
import numpy

#Extension modules

from numpy import *
import cart

class cartModule:

    def __init__(self):
	self.modulename='cart'
        self.cart=cart.cartinterface
	os.chdir('inputs/'+self.modulename)
	self.cart.cart_mpi_init()
	os.chdir('../..')

    def sifInitialize(self,properties,conditions):
	pass
    
    def paramInput(self):
	os.chdir('inputs/'+self.modulename)
	self.cart.cart_param_input()
	os.chdir('../../')

    def initData(self):
	os.chdir('inputs/'+self.modulename)
	self.cart.cart_init_data()
        self.data={'q-data':self.cart.q,
                   's-data':self.cart.rhs,
                   'vol-data':self.cart.vol,
                   'myid-data':self.cart.myid,
                   'n-data':self.cart.ninstances,
                   'jmax-data':self.cart.jmax,
                   'kmax-data':self.cart.kmax,
                   'lmax-data':self.cart.lmax,
                   'h-data':self.cart.h,
                   'freq-data':self.cart.freq,
                   'tcomp_ts-data':self.cart.tcomp_ts,
                   'tcomm_ts-data':self.cart.tcomm_ts,
                   'timecomm-data':self.cart.timecomm,
                   'timerank-data':self.cart.myid_temporal}
	os.chdir('../..')

    def rhs(self,viscous=True,bdf=True):
        self.cart.cart_rhs_inviscid()
        if viscous:
            self.cart.cart_rhs_viscous()
	if bdf:
	    self.cart.cart_bdf_source()

    def lhs(self,it):
        self.cart.cart_lhs(it)

    def update_time(self):
	self.cart.cart_update_time()
        
    def inputParam(self,cfl,bctyp,jmax,kmax,lmax,irhs,ilhs):
        self.cart.cart_init_param(cfl,bctyp,jmax,kmax,lmax,irhs,ilhs)
    
    def computeNorm(self,rres_two,qres_two,rres_inf,qres_inf):
	self.cart.ts_compute_norm(rres_two,qres_two,rres_inf,qres_inf)

    def outputData(self):
	os.chdir('outputs/'+self.modulename)
	self.cart.cart_output()
	os.chdir('../..')

    def finish(self):
	os.chdir('outputs/'+self.modulename)
	self.cart.cart_cleanup()
	os.chdir('../..')
