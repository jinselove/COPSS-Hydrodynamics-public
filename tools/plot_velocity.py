# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 17:45:23 2015

@author: xzhao
"""


import numpy as np
import matplotlib.pyplot as plt
import math

## -------------------------------------------------------------------
## Approximation for the correction factor of Stokes' law
#def correction_factor_bohlin(r_r0):
#    val = 1.0 - 2.10444*r_r0 + 2.08877*math.pow(r_r0,3) \
#          - 0.94813*math.pow(r_r0,5) - 1.372*math.pow(r_r0,6) \
#          + 3.87*math.pow(r_r0,8) - 4.19*math.pow(r_r0,10);
#    return 1.0/val;
#
#
#def correction_factor_haberman(r_r0):
#    f1 = 1.0 - 0.75857*math.pow(r_r0,5);
#    f2 = 1.0 - 2.1050*r_r0 + 2.0865*math.pow(r_r0,3) \
#       - 1.7068*math.pow(r_r0,5) + 0.72603*math.pow(r_r0,6);
#    return f1/f2;
#
#def correction_factor_linear(r_r0):
#    f1 = 3.0*math.pow(r_r0, 10)
#    f2 = 3.0*math.pow(r_r0, 9)
#    return f1/f2;
#    
#
#
## Legend and labels
#vtype = np.array([ "Bohlin's approximate solution", 
#                   "Haberman's approximate solution", 
#                   "GGEM/IBM solution", "Exact Solution" ])
#labelplot = [ str( vt ) for vt in vtype ];
#
#
## analytical solution   
#ndiv = 308;
#r_r0_bo = np.linspace( 0.0, 0.7, ndiv );
#r_r0_ha = np.linspace( 0.0, 0.85, ndiv );
#f_bo = [ correction_factor_bohlin(ri) for ri in r_r0_bo ];
#f_ha = [ correction_factor_haberman(ri) for ri in r_r0_ha ];
#plt.figure(1,figsize=(10,7), dpi=150)
#plt.plot(r_r0_bo, f_bo, 'r-',label=labelplot[0], linewidth=2)  # semilogy/plot
#plt.plot(r_r0_ha, f_ha, 'b-',label=labelplot[1], linewidth=2)
#
## numerical data
#r_r0_num = [ 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ];
#f_num_avg = [ 1.0, 1.305789, 1.752790, 2.508859, 3.861271, 6.330767, 10.115102, 13.049933, 13.961139 ]; 
#f_num_min = [ 1.0, 1.324881, 1.788028, 2.583644, 4.032718, 6.709658, 10.928364, 22.667697, 56.633128 ];
##plt.plot(r_r0_num, f_num_avg, 'bs')
#plt.plot(r_r0_num, f_num_min, 'ko',label=labelplot[2], linewidth=2)
#
#plt.legend(bbox_to_anchor=(0.10, 0.85), loc=2, borderaxespad=0., fontsize=20)
#plt.xlabel(r'$R/R_{0}$', fontsize=20)
#plt.ylabel(r'$\frac{F}{6 \pi \mu R U_{0}}$', fontsize=25)
#plt.show()
#
#plt.savefig('ggem_ibm_validation.eps', format='eps', dpi=2800)



# -------------------------------------------------------------------
## read and plot the history data
#tmpA = np.loadtxt("output_velocity_profile_x_3P.txt");
#mA = tmpA.shape[0]
#nA = tmpA.shape[1]
#
## coordinates of cross section
#x  = tmpA[:,0];      # coordinates
#vx = np.zeros( (mA,4), np.float ); # total | global | local | exact
#vy = np.zeros( (mA,4), np.float ); #           |        |
#vz = np.zeros( (mA,4), np.float ); #          FEM      Green
#
#for i in range(0,4):
#    vx[:,i] = tmpA[:,i*3+1]
#    vy[:,i] = tmpA[:,i*3+2]
#    vz[:,i] = tmpA[:,i*3+3]
#    
#
#colorplot = ["red", "blue", "green", "cyan", "black" ]
#vtype = np.array([ "Total velocity", "Global velocity", "Local velocity", "Exact Solution" ])
#labelplot = [ str( vt ) for vt in vtype ];
#
#
#plt.figure(3,figsize=(10,7), dpi=150)
#for i in range(0, 4):
#    if i<3:
#        plt.plot(x, vx[:,i], label=labelplot[i], color=colorplot[i], linewidth=2);
#    else :
#        plt.plot(x, vx[:,i], label=labelplot[i], color=colorplot[i], linewidth=2, ls='--');        
#        
#plt.xlabel("x", fontsize=25)
#plt.ylabel("u", fontsize=25)
##plt.title("velocity profile in the z-direction")
#plt.legend(bbox_to_anchor=(0.715, 1.0), loc=2, borderaxespad=0.)
#plt.show()

#plt.savefig('ggem_pointforce_validation.eps', format='eps', dpi=1200)


# ---------------------------------------------------------------------
# plot a bounding box
#plt.axvline(x=-1., ymin=0., ymax = 10, linewidth=2, color='k')
#plt.axvline(x=+1., ymin=0., ymax = 10, linewidth=2, color='k')


#plt.figure(1, dpi=100)
#plt.figure(2,figsize=(10,5), dpi=100)
#plt.plot(x,vy, color=colorplot[0], linewidth=2);
#plt.xlabel("x")
#plt.ylabel("v")
#plt.title("velocity profile in the x-direction")
#
#plt.figure(3,figsize=(10,5), dpi=100)
#plt.plot(x,vz, color=colorplot[0], linewidth=2);
#plt.xlabel("x")
#plt.ylabel("w")
#plt.title("velocity profile in the x-direction")
# ---------------------------------------------------------------------



colorplot = ["red", "blue", "green", "cyan", "black" ]
vtype = np.array([ "Practical scaling", "Ideal scaling" ])
labelplot = [ str( vt ) for vt in vtype ];



# -------------------------------------------------------------------
# scalability test:
# domain size (x,y,z)=[-10,+10], mesh 60X60X60, 2.9M dofs totally. 
# rtol=1E-6
# -------------------------------------------------------------------
# msd queue
#ncores = [ 1,      2,     4,     8,     16,    32,    64,   96,   128,  160, 200, 240 ];
#wtime  = [ 1615.7, 931.2, 519.9, 289.4, 254.7, 132.0, 69.4, 47.8, 37.3, 30.6,26.1,22.4 ];

# regular node of blues(8, 658.8/706.0) (9, 503.8) (320, 14.0)
ncores = [ 1,      2,     4,     16,    32,    64,   128,  256, 320, 480 ];
wtime  = [ 1470.2, 989.7, 674.1, 191.5, 99.9, 52.7, 28.2, 16.2, 14.0, 10.2 ];


itime  = np.zeros( len(ncores) );
for i in range(0, len(ncores)):
    itime[i] = wtime[0]/ncores[i];


# --- Plot scalability ---semilogx, semilogy, loglog, plot
plt.figure(1,figsize=(8,5), dpi=150)
plt.loglog(ncores,wtime, 'rs-', label=labelplot[0],linewidth=2);
plt.loglog(ncores,itime, 'r--', label=labelplot[1],linewidth=2);
plt.xlabel("Number of CPUs")
plt.ylabel("Computational time (s)")
plt.legend()
plt.show()
plt.savefig('scalability.eps', format='eps', dpi=150)


# -------------------------------------------------------------------
# speedup of N points, run on 4xppn=16 =>64 cores
# -------------------------------------------------------------------
#N      = [  100,   1000,  2000,  4000,  6000,  8000,  10000, 12500,  15000 ];
#
## # # alpha*h = 0.6
#Ntime1 = [  113.7, 164.1, 185.6, 340.1, 573.9, 889.5, 1335.0,1962.35,2807.9 ];
#Ntime2 = [  105.6, 118.6, 136.7, 211.6, 325.9, 478.2, 716.1, 1005.9, 1435.7 ];
##   iter    17,    16,    16     16     16     16     16     16     16  
#
## # # alpha*h = 1.0
##Ntime1 = [  116.8, 174.5, 188.9, 339.4, 569.1, 895.1, 1319.5,1986.1,2837.9 ];
##Ntime2 = [  108.8, 136.1, 143.1, 218.0, 324.2, 497.1, 689.8, 978.9, 1423.4 ];
###   iter    17,    17,    17     17     16     16     16     16     16  
#
#ntime  = np.zeros( len(N) );
#for i in range(0, len(N)):
##    ntime[i] = Ntime2[0]*math.exp( (1.8E-4)*N[i] );
#    ntime[i] = 0.01*N[i]*math.log( N[i] );
#    
#
## --- Plot N-dependence ---
#vtype = np.array([ "The 1st solve", "The 2nd solve" ])
#labelplot = [ str( vt ) for vt in vtype ];
#plt.figure(2,figsize=(10,4), dpi=100)
##plt.plot(N,Ntime1, 'bo-', label=labelplot[0],linewidth=2);
#plt.plot(N,Ntime2, 'ks-', label=labelplot[1],linewidth=2);
#plt.plot(N,ntime, 'k--', label=labelplot[1],linewidth=2);
#plt.xlabel("Number of points (N)")
#plt.ylabel("Computational time (s)")
#plt.legend(bbox_to_anchor=(0.5, 1.0), loc=2, borderaxespad=0.)
#plt.show()





# -------------------------------------------------------------------
# convergence test:
# domain size 10X10X10 
# -------------------------------------------------------------------
#alpha = [ 0.05, 0.10,  0.20, 0.30, 0.40, 0.5,
#            1., 2., 3., 4];
#L2_norm = [ 1.476035E-10,3.746274E-09,5.802349E-08,2.117199E-07,6.224134E-07,
#1.350082E-06,4.271880E-04,1.760053E-02,3.142662E-02,3.864256E-02
#];
#
#
#alpha_1 = [ 0.05, 0.10, 0.15, 0.20, 0.25, 
#            0.30, 0.35, 0.40, 0.45, 0.5,
#            1., 2., 3., 4., 5., 6., 7., 8.];
#
## mesh size 10X10X10
#L2_norm_1 =[ 3.565112E-08,8.142992E-07,4.317537E-06,1.212812E-05,
#2.497848E-05,4.620671E-05,8.435798E-05,1.434607E-04,
#2.239432E-04,3.222335E-04,2.943570E-03,2.213465E-02,
#5.017971E-02,1.614902E-01,2.320975E-01,1.749748E-01,
#4.324703E-02,1.451656E-01 ];
#
## mesh size 20X20X20
#L2_norm_2 = [2.761795E-09,5.585424E-08,3.061887E-07,8.614827E-07,
#1.746009E-06,3.183387E-06,5.679496E-06,9.503402E-06,
#1.467817E-05,2.092603E-05,5.212515E-04,1.760893E-02,
#3.123014E-02,3.760233E-02,3.999702E-02,6.288481E-02,
#1.160841E-01,1.876017E-01 ]; 
#    
#
#mesh_size = np.array([ "10X10X10", "20X20X20", "40X40X40", "80X80X80" ])
#label_mesh = [ str( ms ) for ms in mesh_size ];
#
#plt.figure(1001)
#plt.loglog(alpha_1[0:12],L2_norm_1[0:12],'ro-',label=label_mesh[0],linewidth=2)
#plt.loglog(alpha_1[0:12],L2_norm_2[0:12],'b^-',label=label_mesh[1],linewidth=2)
#plt.loglog(alpha,L2_norm,'ks-',label=label_mesh[2],linewidth=2)
#plt.xlabel(r'$\alpha$')
#plt.ylabel("L1 error")
#plt.legend()
#plt.show()
#
## plot the slope
#x = [ 0.1, 1. ]
#y1 = [1E-9, 1E-8]
#y2 = [1E-9, 1E-7]
#y3 = [1E-9, 1E-6]
#plt.loglog(x,y1,'g.-',linewidth=1)
#plt.loglog(x,y2,'g.-',linewidth=1)
#plt.loglog(x,y3,'g.-',linewidth=1)
