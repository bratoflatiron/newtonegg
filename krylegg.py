import numpy as np 
import tkinter
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from scipy import optimize
from scipy.integrate import odeint 
from scipy.integrate import solve_ivp 
from scipy.io import savemat
import newkry
#print (newkry.__doc__)

# #--------------------------ALL THIS IS PREPERATION-------------------------------------
# def prepare():
#     ## Load the basic constants of the problem 
#     [a1,a2,a3,rad,delta, nquad, npols, finout] = newkry.loadconstant()

#     ## Define the constants for the surfaces
#     rmax = rad*0.6;  rmSL = (rad-delta)*0.6; 
#     a1SL = a1-delta; a2SL = a2-delta; a3SL = a3-delta; # Aspect ratio for the slip surface
#     ifc  = int(0);

#     nptchB         = newkry.get_ellipsoid_mem(a1,   a2,   a3,   rmax, ifc) 
#     [rmsl,nptchSL] = newkry.rmxestimate(nptchB, ifc, a1SL, a2SL, a3SL)
#     nptB = nptchB*npols; nptSL = nptchSL*npols;

#     ## Compute surface geometries and quadrature weights
#     [nordB, ixyzB, iptB, surfB,   srBcoeff,   wtsB]   = newkry.boundb (rmax,nquad,nptchB, nptB)
#     [nordSL,ixyzSL,iptSL,surfslip,srslipcoeff,wtsslip]= newkry.boundsl(rmSL,nquad,nptchSL,nptSL)

#     ## Obtain points for interpolation and matrix for interpolation
#     [xinterp,matinterp] = newkry.get_interp_mat(surfslip,rmax,nptB,npols,nordB,ixyzB,iptB,
#                                             nptSL,nptchB)

#     ## Load precomputed data 
#     pdir = "/mnt/home/bchakrabarti/Documents/oocyte/BIE_V2/pdataell/a10_10/"
#     [LB_grad,tracD_mat] = newkry.readprecomp(nptB,nptSL,pdir)

#     ## Define the Jmat for gmres. Jmat = LB_grad x matinterp
#     Jmat = np.matmul(LB_grad,matinterp)

#     return nptB, nptSL, nptchB, nordB,  ixyzB, iptB, rmax, surfB, xinterp, LB_grad, Jmat 
#--------------------------------------------------------------------------------------
#--------------------------INITIAL CONDITIONS AND DENSITY------------------------------

# ## Load data/precomputation
# [nptB, nptSL, nptchB, nordB,  ixyzB, iptB, rmax, surfB, xinterp, LB_grad, Jmat] = prepare()

# ## Obtain initial condition 
# choice = int(0);
# [nDB,nMT] = newkry.restart(choice,surfB,nptB)

# ## Obtain surface density
# rho = newkry.density(surfB,'no',nptB)
# #--------------------------TIME MARCHING------------------------------
# # brink = int(1); delt = 5e-3; steps = int(400);
# # sol = newkry.timemarch(delt,rmax,nDB,rho,xinterp,LB_grad,Jmat,nMT, 
# #                 brink,steps,nordB, ixyzB, iptB,nptB, nptSL,nptchB)
# #----------------------RESIDUE FUNCTION----------------------------
# brink = int(1); delt = 5e-3; steps = int(5);
# def residue(u): 
#     init = np.reshape(u,(3,nptB)); 
#     sol = newkry.timemarch(delt,rmax,nDB,rho,xinterp,LB_grad,Jmat,init, 
#                  brink,steps,nordB, ixyzB, iptB,nptB, nptSL,nptchB)

#     sol = np.reshape(sol,(3*nptB, 1)); 
#     res = u-sol;
#     return res

# guess = np.reshape(nMT,(3*nptB, 1)); 
# err = residue(guess); print("initial error:", np.linalg.norm(err))
# print("-------starting newton-krylov--------")
# sol = optimize.newton_krylov(residue, guess, f_tol = 1e-5, x_tol = 1e-4, verbose=True); 
# sol = np.reshape(sol,(3,nptB))
# savemat("krylovsol.mat", mdict={'sol': sol})