import numpy as np
import math
import json
from collections import Counter
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import math
import os, sys
import datetime
from sys import platform
from sys import stdout
from global_param import *
from global_param2 import *
from global_func import *
from read_structure import *
from user_input import *

#~~~~~~~~~~~~~~~~~~~~~~~~~Create an exclusion list~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

exclusion_list = []


#~~~~~~~~~~~~~~~~~~~~~~~Bonded Interactions within DNA~~~~~~~~~~~~~~~~~~~~~#

bond_dna = []      # Format------[[index1, index2, k_spring, dist_equilibrium],[],[]]
for j in range (0,ndna):
    for i in range (0,dict_nres["nres_tot_chain_A"]):
        S_bead_index = 3*i + j*dict_nbead["bead_tot_chain_A"]
        P_next_bead_index = S_bead_index + 2 
        P_prev_bead_index = S_bead_index - 1
        B_bead_index = S_bead_index + 1 

        if P_next_bead_index in range(0 + j*dict_nbead["bead_tot_chain_A"], dict_nbead["bead_tot_chain_A"] + j*dict_nbead["bead_tot_chain_A"]):
            key1=str(atom_type[S_bead_index])+str(atom_type[P_next_bead_index])
            bond_dna.append(tuple([S_bead_index, P_next_bead_index, bond_type[key1], calculate_dist(S_bead_index, P_next_bead_index)]))
            exclusion_list.append(tuple([S_bead_index, P_next_bead_index]))
        if P_prev_bead_index in range(0 + j*dict_nbead["bead_tot_chain_A"],dict_nbead["bead_tot_chain_A"] + j*dict_nbead["bead_tot_chain_A"]):
            key2=str(atom_type[P_prev_bead_index])+str(atom_type[S_bead_index])
            bond_dna.append(tuple([S_bead_index, P_prev_bead_index, bond_type[key2], calculate_dist(P_prev_bead_index, S_bead_index)]))
            exclusion_list.append(tuple([S_bead_index, P_prev_bead_index]))
        key3=str(atom_type[S_bead_index])+str(residue[B_bead_index])    
        bond_dna.append(tuple([S_bead_index, B_bead_index, bond_type[key3], calculate_dist(S_bead_index, B_bead_index)]))
        exclusion_list.append(tuple([S_bead_index, B_bead_index]))
        
print(len(bond_dna))



#~~~~~~~~~~~~~~~~~~~~~~~Angle within DNA~~~~~~~~~~~~~~~~~~~~~~~~~~#


angle_dna = []   # Format------[[index1, index2, index3, k_spring, angle_equilibrium],[],[]]


for j in range (0,ndna):
    for i in range (0, dict_nres["nres_tot_chain_A"]):
        S_bead_index = 3*i + j*dict_nbead["bead_tot_chain_A"]
        P_next_bead_index = 3*i + 2 + j*dict_nbead["bead_tot_chain_A"]
        P_prev_bead_index = 3*i - 1 + j*dict_nbead["bead_tot_chain_A"]
        B_bead_index = 3*i + 1 + j*dict_nbead["bead_tot_chain_A"]

        if P_prev_bead_index in range(0 + j*dict_nbead["bead_tot_chain_A"],dict_nbead["bead_tot_chain_A"] + j*dict_nbead["bead_tot_chain_A"]):
            key1=str(atom_type[P_prev_bead_index]) + str(atom_type[S_bead_index]) + str(residue[B_bead_index])
            angle_dna.append(tuple([P_prev_bead_index, S_bead_index, B_bead_index, angle_type[key1], calculate_angle(P_prev_bead_index, S_bead_index, B_bead_index)]))
            exclusion_list.append(tuple([P_prev_bead_index, B_bead_index]))
        if P_next_bead_index in range(0 + j*dict_nbead["bead_tot_chain_A"],dict_nbead["bead_tot_chain_A"] + j*dict_nbead["bead_tot_chain_A"]):    
            key2 = str(str(residue[B_bead_index]) + str(atom_type[S_bead_index]) + str(atom_type[P_next_bead_index]))
            angle_dna.append(tuple([B_bead_index, S_bead_index, P_next_bead_index, angle_type[key2], calculate_angle(B_bead_index, S_bead_index, P_next_bead_index)]))
            exclusion_list.append(tuple([B_bead_index, P_next_bead_index]))

        if P_next_bead_index in range(0 + j*dict_nbead["bead_tot_chain_A"],dict_nbead["bead_tot_chain_A"] + j*dict_nbead["bead_tot_chain_A"]) and P_prev_bead_index in range(0 + j*dict_nbead["bead_tot_chain_A"],dict_nbead["bead_tot_chain_A"] + j*dict_nbead["bead_tot_chain_A"]):
            key3 = str(atom_type[P_prev_bead_index]) + str(atom_type[S_bead_index]) + str(atom_type[P_next_bead_index])
            angle_dna.append(tuple([P_prev_bead_index, S_bead_index, P_next_bead_index, angle_type[key3], calculate_angle(P_prev_bead_index, S_bead_index, P_next_bead_index)]))
            exclusion_list.append(tuple([P_prev_bead_index, P_next_bead_index]))

print(len(angle_dna))    



#~~~~~~~~~~~~~~~~~~~~Adjacent stacking~~~~~~~~~~~~~~~~~~~~~~#


stacked = []   #Format: S1, P1, P2, S2, P3, B1, B2, U^0_S, l_eql, phi1_eql, phi2_eql #

for j in range (0,ndna):
    for i in range (1, dict_nres["nres_tot_chain_A"]-2):
        S_bead_index = 3*i + j*dict_nbead["bead_tot_chain_A"]
        P_bead_index = S_bead_index - 1
        P_next_bead_index = P_bead_index + 3
        S_next_bead_index = S_bead_index + 3

        P_next_next_bead_index = P_next_bead_index + 3

        B_bead_index = S_bead_index + 1

        B_next_bead_index = B_bead_index + 3

        key = str(residue[B_bead_index]) + str(residue[B_next_bead_index])

        U_0 = -float(stack_type[key][0]) + k_B*(temp - float(stack_type[key][2]))*float(stack_type[key][1])  #U_0 = -h + k+{B}*(T - T_m)*s

        stacked.append(tuple([S_bead_index, P_bead_index, P_next_bead_index, S_next_bead_index, P_next_next_bead_index, B_bead_index, B_next_bead_index, U_0, calculate_dist(B_bead_index, B_next_bead_index), calculate_dihedral(P_bead_index, S_bead_index, P_next_bead_index, S_next_bead_index), calculate_dihedral(P_next_next_bead_index, S_next_bead_index, P_next_bead_index, S_bead_index)]))

print(len(stacked))



#~~~~~~~~~~~~~~~~~~~WC Hydrogen bonding~~~~~~~~~~~~~~~~~~~~~~~~~~#


hbonded = []     #~~~~~~~~~~~~Format: B1, S1, P2, B5, S5, P6, r, theta1, theta2, phi1, phi2, phi3, U_0~~~~~~~~~~~~~#
sum_tot = 0.0
for i in range (1, dict_nres["nres_tot_chain_A"]-1):
    B_bead_index_st1 = 3*i + 1
    B_bead_index_st2 = 3*(dict_nres["nres_tot_chain_A"] + dict_nres["nres_tot_chain_B"] -1 - i)

    S_bead_index_st1 = B_bead_index_st1 - 1
    P_bead_index_st1 = S_bead_index_st1 + 2

    S_bead_index_st2 = B_bead_index_st2 -1
    P_bead_index_st2 = S_bead_index_st2 + 2

    key = str(residue[B_bead_index_st1]) + str(residue[B_bead_index_st2])
    if key in hbond_type:
        hbonded.append(tuple([B_bead_index_st1, S_bead_index_st1, P_bead_index_st1, B_bead_index_st2, S_bead_index_st2, P_bead_index_st2, calculate_dist(B_bead_index_st1, B_bead_index_st2), calculate_angle(S_bead_index_st1, B_bead_index_st1, B_bead_index_st2), calculate_angle(S_bead_index_st2, B_bead_index_st2, B_bead_index_st1), calculate_dihedral(S_bead_index_st1, B_bead_index_st1, B_bead_index_st2, S_bead_index_st2), calculate_dihedral(P_bead_index_st2, S_bead_index_st2, B_bead_index_st2, B_bead_index_st1), calculate_dihedral(P_bead_index_st1, S_bead_index_st1, B_bead_index_st1, B_bead_index_st2), hbond_type[key]]))

    else:
        print("key {key} not found in Hbond list\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~Fene bond in SOP-SC model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

bond_pro = []        #~~~~~~~~~~~~~Format: index1, index2, r_cry~~~~~~~~~~~~~~~~~#

for i in range (2, nchain):
   tot_residue = int(info_dna_pro[name_list[i]]/2)
   for j in range (0, tot_residue):
       bb_ind = get_pro_bb_index(i, j)
       sc_ind = get_pro_sc_index(i, j)
       bond_pro.append(tuple([bb_ind, sc_ind, calculate_dist(bb_ind, sc_ind)]))
       exclusion_list.append(tuple([bb_ind, sc_ind]))
       if j < (tot_residue - 1): 
          bb_ind_nxt = get_pro_bb_index(i, j+1)
          bond_pro.append(tuple([bb_ind, bb_ind_nxt, calculate_dist(bb_ind, bb_ind_nxt)]))
          exclusion_list.append(tuple([bb_ind, bb_ind_nxt]))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Angle in SOP-SC model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

angle_pro = []       #~~~~~~~~~~~~~~Format: index1, index2, sig(index1) + sig(index2)~~~~~~~~~~~~~~#

for i in range (2, nchain):
   tot_residue = int(info_dna_pro[name_list[i]]/2)
   for j in range (0, tot_residue):
       bb_ind = get_pro_bb_index(i, j)
       sc_ind = get_pro_sc_index(i, j)
       if j < (tot_residue - 1):
           bb_ind_nxt = get_pro_bb_index(i, j+1)
           sc_ind_nxt = get_pro_sc_index(i, j+1)
           angle_pro.append(tuple([bb_ind, sc_ind_nxt, radii_scale*(pro_bb_radius + get_radii_pro(i, j+1))]))
           angle_pro.append(tuple([sc_ind, bb_ind_nxt, radii_scale*(get_radii_pro(i, j) + pro_bb_radius)]))
           exclusion_list.append(tuple([bb_ind, sc_ind_nxt]))
           exclusion_list.append(tuple([sc_ind, bb_ind_nxt]))
           if j < (tot_residue - 2):
               bb_ind_nxt_nxt = get_pro_bb_index(i, j+2)
               angle_pro.append(tuple([bb_ind, bb_ind_nxt_nxt, 2*pro_bb_radius]))
               exclusion_list.append(tuple([bb_ind, bb_ind_nxt_nxt]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~List of Native Interactions within proteins~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

native_pro = []     #~~~~~~~~~~~~~~~~Format: particle1, particle2, r12, eps12~~~~~~~~~~~~~~~~~~~~~~~~#


for i in range (2, nchain):
   tot_residue = int(info_dna_pro[name_list[i]]/2)
   for j in range (0, tot_residue - 1):
      bb_ind_j = get_pro_bb_index(i, j)
      for k in range (j+1, tot_residue):
         bb_ind_k = get_pro_bb_index(i, k)
         d = get_dist_along_chain(bb_ind_j, bb_ind_k)
         if d >= bond_sep_nat:
            dist = calculate_dist(bb_ind_j, bb_ind_k)
            if dist <= pro_nat_cut:
                native_pro.append(tuple([bb_ind_j, bb_ind_k, dist, pro_eps_bb]))
                exclusion_list.append(tuple([bb_ind_j, bb_ind_k]))
for i in range (2, nchain):
   tot_residue = int(info_dna_pro[name_list[i]]/2)
   for j in range (0, tot_residue - 1):
      bb_ind_j = get_pro_bb_index(i, j)
      for k in range (j+1, tot_residue):
         sc_ind_k = get_pro_sc_index(i, k)
         d = get_dist_along_chain(bb_ind_j, sc_ind_k)
         if d >= bond_sep_nat:
            dist = calculate_dist(bb_ind_j, sc_ind_k)
            if dist <= pro_nat_cut:
                native_pro.append(tuple([bb_ind_j, sc_ind_k, dist, pro_eps_bs]))
                exclusion_list.append(tuple([bb_ind_j, sc_ind_k]))


for i in range (2, nchain):
   tot_residue = int(info_dna_pro[name_list[i]]/2)
   for j in range (0, tot_residue - 1):
      sc_ind_j = get_pro_sc_index(i, j)
      for k in range (j+1, tot_residue):
         bb_ind_k = get_pro_bb_index(i, k)
         d = get_dist_along_chain(sc_ind_j, bb_ind_k)
         if d >= bond_sep_nat:
            dist = calculate_dist(sc_ind_j, bb_ind_k)
            if dist <= pro_nat_cut:
               native_pro.append(tuple([sc_ind_j, bb_ind_k, dist, pro_eps_bs]))
               exclusion_list.append(tuple([sc_ind_j, bb_ind_k]))

for i in range (2, nchain):
   tot_residue = int(info_dna_pro[name_list[i]]/2)
   for j in range (0, tot_residue - 1):
      sc_ind_j = get_pro_sc_index(i, j)
      for k in range (j+1, tot_residue):
         sc_ind_k = get_pro_sc_index(i, k)
         d = get_dist_along_chain(sc_ind_j, sc_ind_k)
         if d == bond_sep_nat:
            exclusion_list.append(tuple([sc_ind_j, sc_ind_k]))   # Adjacent side-chains do not see each other#
         if d >= bond_sep_nat_sc:
            dist = calculate_dist(sc_ind_j, sc_ind_k)
            if dist <= pro_nat_cut:

               residue1 = residue[sc_ind_j]
               residue2 = residue[sc_ind_k]

               eps_SS = pro_eps_ss*abs(0.7 - get_eps_from_BT_pot(residue1,residue2))
               native_pro.append(tuple([sc_ind_j, sc_ind_k, dist, eps_SS]))
               exclusion_list.append(tuple([sc_ind_j, sc_ind_k]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~List of native interactions between different chains~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#BB-BB
for i in range (2, nchain):
    nres_chain_i = int(info_dna_pro[name_list[i]]/2)
    for i_res in range (0, nres_chain_i):
        bb_ind_i = get_pro_bb_index(i, i_res)
        for j in range (i+1, nchain):
            nres_chain_j = int(info_dna_pro[name_list[j]]/2)
            for j_res in range (0, nres_chain_j):
                bb_ind_j = get_pro_bb_index(j, j_res)
                dist = calculate_dist(bb_ind_i, bb_ind_j)
                if dist <= pro_nat_cut:
                    native_pro.append(tuple([bb_ind_i, bb_ind_j, dist, pro_eps_bb]))
                    exclusion_list.append(tuple([bb_ind_i, bb_ind_j]))
                     
#BB-SC
for i in range (2, nchain):
    nres_chain_i = int(info_dna_pro[name_list[i]]/2)
    for i_res in range (0, nres_chain_i):
        bb_ind_i = get_pro_bb_index(i, i_res)
        for j in range (i+1, nchain):
            nres_chain_j = int(info_dna_pro[name_list[j]]/2)
            for j_res in range (0, nres_chain_j):
                sc_ind_j = get_pro_sc_index(j, j_res)
                dist = calculate_dist(bb_ind_i, sc_ind_j)
                if dist <= pro_nat_cut:
                    native_pro.append(tuple([bb_ind_i, sc_ind_j, dist, pro_eps_bs]))
                    exclusion_list.append(tuple([bb_ind_i, sc_ind_j]))

for i in range (2, nchain):
    nres_chain_i = int(info_dna_pro[name_list[i]]/2)
    for i_res in range (0, nres_chain_i):
        sc_ind_i = get_pro_sc_index(i, i_res)
        for j in range (i+1, nchain):
            nres_chain_j = int(info_dna_pro[name_list[j]]/2)
            for j_res in range (0, nres_chain_j):
                bb_ind_j = get_pro_bb_index(j, j_res)
                dist = calculate_dist(sc_ind_i, bb_ind_j)
                if dist <= pro_nat_cut:
                    native_pro.append(tuple([sc_ind_i, bb_ind_j, dist, pro_eps_bs]))
                    exclusion_list.append(tuple([sc_ind_i, bb_ind_j]))

#SC-SC
for i in range (2, nchain):
    nres_chain_i = int(info_dna_pro[name_list[i]]/2)
    for i_res in range (0, nres_chain_i):
        sc_ind_i = get_pro_sc_index(i, i_res)
        for j in range (i+1, nchain):
            nres_chain_j = int(info_dna_pro[name_list[j]]/2)
            for j_res in range (0, nres_chain_j):
                sc_ind_j = get_pro_sc_index(j, j_res)
                dist = calculate_dist(sc_ind_i, sc_ind_j)
                if dist <= pro_nat_cut:
                    residue1 = residue[sc_ind_i]
                    residue2 = residue[sc_ind_j]
                    eps_SS = pro_eps_ss*abs(0.7 - get_eps_from_BT_pot(residue1,residue2))
                    native_pro.append(tuple([sc_ind_i, sc_ind_j, dist, eps_SS]))
                    exclusion_list.append(tuple([sc_ind_i, sc_ind_j]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DNA-PRO Native interactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dna_pro_native = []         #~~~~~~~~~~~~~~~~~~~~Format: particle1, particle2, r12, eps12~~~~~~~~~~~~~~~~~~~~~~~~#
                            #NOTE: P bead of DNA does not interact with SC beads of proteins except for electrostatics

for i in range (0, bead_tot_dna):
   for j in range (bead_tot_dna, bead_tot):
      dist = calculate_dist(i, j)
      if dist < pro_dna_nat_cut:
         if atom_type[i] == "DP2" and atom_type[j] == "BK":
            eps = scale_nat_dna_pro*P_BB[residue[j]][residue[i]]
            dna_pro_native.append(tuple([i, j, dist, eps])) 
            exclusion_list.append(tuple([i,j]))
         if atom_type[i] == "DS2" and atom_type[j] == "BK":
            eps = scale_nat_dna_pro*S_BB[residue[j]][residue[i]]
            dna_pro_native.append(tuple([i, j, dist, eps]))
            exclusion_list.append(tuple([i,j]))
         if atom_type[i] == "DS2" and atom_type[j] == "SC":
            eps = scale_nat_dna_pro*S_SC[residue[j]][residue[i]]
            dna_pro_native.append(tuple([i, j, dist, eps]))
            exclusion_list.append(tuple([i,j]))
         if atom_type[i] == "DB2" and atom_type[j] == "BK":
            eps = scale_nat_dna_pro*B_BB[residue[j]][residue[i]]
            dna_pro_native.append(tuple([i, j, dist, eps]))
            exclusion_list.append(tuple([i,j]))
         if atom_type[i] == "DB2" and atom_type[j] == "SC":
            eps = scale_nat_dna_pro*B_SC[residue[j]][residue[i]]
            dna_pro_native.append(tuple([i, j, dist, eps]))
            exclusion_list.append(tuple([i,j]))

#~~~~~~~~~~~~~~~~~~~~~~~~~DNA-DNA Excluded Volume Interactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dna_dna_ex_vol = []

for i in range (0, bead_tot):
    if i < bead_tot_dna:
        eps = 1.0
        dna_dna_ex_vol.append(tuple([eps]))

    else:
        eps = 0.0
        dna_dna_ex_vol.append(tuple([eps]))

#~~~~~~~~~~~~~~~~~~~~~~~~~PRO-PRO Excluded Volume Interactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pro_pro_ex_vol = []

for i in range (0, bead_tot):
    if i < bead_tot_dna:
       eps = 0.0
       rad = 0.0
       pro_pro_ex_vol.append(tuple([eps, rad]))
    else:
        eps = 1.0
        if atom_type[i] == "BK":
            rad = pro_bb_radius
            pro_pro_ex_vol.append(tuple([eps, rad]))
        else:
            rad = radii_pro[residue[i]]
            pro_pro_ex_vol.append(tuple([eps, rad]))

#~~~~~~~~~~~~~~~~~~~~~~~DNA-PRO Excluded Volume Interactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#NOTE: P bead of DNA does not interact with SC beads of proteins except for electrostatics and excluded volume#
#Exclusion list to exclude all P and SC Interactions#


#Between DNA and PRO:

dna_pro_ex_vol = []

for i in range (0, bead_tot):
   if i < bead_tot_dna:
      eps = 1.0
      scale = 1
      if atom_type[i] == "DP2" or atom_type[i] == "DS2":
          dna_pro_ex_vol.append(tuple([dna_bead[atom_type[i]][1], eps, scale]))
      else:
          dna_pro_ex_vol.append(tuple([dna_bead[residue[i]][1], eps, scale]))
   else:
      eps = 1.0
      scale = -1
      if atom_type[i] == "BK":
         dna_pro_ex_vol.append(tuple([pro_bb_radius, eps, scale]))
      else:
         dna_pro_ex_vol.append(tuple([radii_pro[residue[i]], eps, scale])) 






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Electrostatics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#NOTE: Electrostatic interactions within DNA considers counter ion condensation and takes a reduced charge on P (-0.6)
dna_dna_DH = []

for i in range (0, bead_tot):
    if i < bead_tot_dna:
       if atom_type[i] == "DP2":
          dna_dna_DH.append(tuple([norm_charge_p]))
       else:
          dna_dna_DH.append(tuple([0.0]))
    else:
       dna_dna_DH.append(tuple([0.0]))
 
#NOTE: Charges on Protein side-chains are taken to be 1 or -1 depending on the residue
pro_pro_DH = []

for i in range (0, bead_tot):
    if i < bead_tot_dna:
       pro_pro_DH.append(tuple([0.0]))
    else:
       if atom_type[i] == "SC":
          if residue[i] == "LYS":
             pro_pro_DH.append(tuple([1.0]))
          elif residue[i] == "ARG":
             pro_pro_DH.append(tuple([1.0]))
          elif residue[i] == "ASP":
             pro_pro_DH.append(tuple([-1.0]))
          elif residue[i] == "GLU":
             pro_pro_DH.append(tuple([-1.0]))
          else:
             pro_pro_DH.append(tuple([0.0]))
       else:
          pro_pro_DH.append(tuple([0.0]))
#NOTE: Charges on DNA phosphate is taken to be -1 and protein charges are taken to be +1 or -1
dna_pro_DH = []

for i in range (0,bead_tot):
    if i < bead_tot_dna:
       scale=1
       if atom_type[i] == "DP2":
          dna_pro_DH.append(tuple([-1.0, scale]))
       else:
          dna_pro_DH.append(tuple([0.0, scale]))
    else:
       scale = -1
       if atom_type[i] == "SC":
          if residue[i] == "LYS":
             dna_pro_DH.append(tuple([1.0, scale]))
          elif residue[i] == "ARG":
             dna_pro_DH.append(tuple([1.0, scale]))
          elif residue[i] == "ASP":
             dna_pro_DH.append(tuple([-1.0, scale]))
          elif residue[i] == "GLU":
             dna_pro_DH.append(tuple([-1.0, scale]))
          else:
             dna_pro_DH.append(tuple([0.0, scale]))
       else:
          dna_pro_DH.append(tuple([0.0, scale]))    

