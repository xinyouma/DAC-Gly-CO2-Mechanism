import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import time
o_debug = False
R_0=1.2 
NN=9 
MM=18 
D_0=0.2

a  = NN
r0 = R_0
d0 = D_0

m_N = 14
m_H = 1

# loading a trajectory
natom = 323 
filename='temp_traj.xyz'
if o_debug: print (" loading trajectory from file :" + filename)
tic = time.perf_counter()
traj=np.loadtxt(filename, usecols=(1,2,3))
toc = time.perf_counter()
if o_debug: print (f" finish reading traj file in {toc - tic:0.4f} seconds")
num_step = int(len(traj)/natom)


if o_debug: print ("in total " + str(num_step) + " steps")
traj_mat=np.reshape(traj, (num_step, natom,3))

# re-defining/re-calculating the CV
#cn0: COORDINATION GROUPA=11 GROUPB=2,3,5,6,12,13,16,25,125-323 

# list of all hydrogen atoms indecies
list_H = np.array([2,3,5,6,12,13,16,25])
list_H = np.r_[list_H, range(125,324)]
list_H = list_H - 1 # convert list from atomic series to python series

# list of all hydrogen atoms
list_N = np.array(11)
list_N = list_N - 1


for i_step in range(num_step):
    #list_atoms=np.r_[list_N,list_H]
    N_xyz=traj_mat[i_step,list_N,:]
    H_xyz=traj_mat[i_step,list_H,:]
    N_xyz=np.reshape(N_xyz,(1,3))
    
    # getting all N-H distances with respect to N in CO2NH2-CH2-COO-
    r_NH_all = distance_matrix(N_xyz, H_xyz)
    r_NH_all_x = H_xyz[:,0] - N_xyz[:,0]
    r_NH_all_y = H_xyz[:,1] - N_xyz[:,1]
    r_NH_all_z = H_xyz[:,2] - N_xyz[:,2]
    
    rij_d0_r0 = ( r_NH_all - d0 )/r0
    dsdr = -a/r0 * np.power(rij_d0_r0, a-1) / np.power(1 + np.power(rij_d0_r0, a), 2)
    dsdr2= np.power(dsdr,2)
    
    drdx = (H_xyz[:,0] - N_xyz[:,0])/rij_d0_r0
    drdy = (H_xyz[:,1] - N_xyz[:,1])/rij_d0_r0
    drdz = (H_xyz[:,2] - N_xyz[:,2])/rij_d0_r0
    
    drdx2 = np.power(drdx,2)
    drdy2 = np.power(drdy,2)
    drdz2 = np.power(drdz,2)
    
    #f'_i : dsdr
    #f'_i^2 : dsdr2
    #np.power(np.sum(dsdr*drdx)
    
    # \LAMBDA_auto 
    ## Zn = ( \LAMBDA_auto + \LAMBDA_cross ) / \mu
    # \LAMBDA_auto = \sum_i^{N-1} f'_i^2
    # \LAMBDA_cross = 2*\mu/m_H * \sum_i^{N-2}\sum_j^{N-1} f'_i*f'_j/(r_i*r_j)*(r_i_x*r_j_x + r_i_y*r_j_y + r_i_z*r_j_z)
    
    Zn_H = np.sum(dsdr2) / m_H
    dsdr_dsdr_cross = np.outer(dsdr,dsdr)
    ri = r_NH_all[0] 
    ri_rj_outer = np.outer(ri, ri)
    ri_rj_outer_x = np.outer(r_NH_all_x, r_NH_all_x)
    ri_rj_outer_y = np.outer(r_NH_all_y, r_NH_all_y)
    ri_rj_outer_z = np.outer(r_NH_all_z, r_NH_all_z)
    Zn_N = np.sum(dsdr_dsdr_cross*( ri_rj_outer_x + ri_rj_outer_y + ri_rj_outer_z )/ri_rj_outer/m_N)
    
    print(Zn_N + Zn_H)



