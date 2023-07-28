import numpy as np
from itertools import groupby
from global_param import *
from global_param2 import *
from read_structure import *
from user_input import *

#Note: Function to read json input files if any 

def read_file(fname):
    with open(fname, "r") as F:
        data = F.read()
        s = json.loads(data)
        return s

#Note: Calculate 3D distance between two points

def calculate_dist(ind1, ind2):
    rx = float(coor_x[ind1]) - float(coor_x[ind2])
    ry = float(coor_y[ind1]) - float(coor_y[ind2])
    rz = float(coor_z[ind1]) - float(coor_z[ind2])
    r = math.sqrt(rx*rx + ry*ry + rz*rz)
    return r

#Note: Calculate angle in radians for a set of three points

def calculate_angle(ind1, ind2, ind3):
    a = np.array([float(coor_x[ind1]), float(coor_y[ind1]),float(coor_z[ind1])])
    b = np.array([float(coor_x[ind2]), float(coor_y[ind2]),float(coor_z[ind2])])
    c = np.array([float(coor_x[ind3]), float(coor_y[ind3]),float(coor_z[ind3])])

    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return angle

#Note: Calculate dihedral in radians for a set of four points

def calculate_dihedral(ind1, ind2, ind3, ind4):
    u1 = np.array([float(coor_x[ind1]), float(coor_y[ind1]),float(coor_z[ind1])])
    u2 = np.array([float(coor_x[ind2]), float(coor_y[ind2]),float(coor_z[ind2])])
    u3 = np.array([float(coor_x[ind3]), float(coor_y[ind3]),float(coor_z[ind3])])
    u4 = np.array([float(coor_x[ind4]), float(coor_y[ind4]),float(coor_z[ind4])])


    a1 = u2 - u1
    a2 = u3 - u2
    a3 = u4 - u3

    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm

    return rad

#Note: Get bead index of SOP-SC sidechain bead

def get_pro_sc_index(serial, residue_id):
   if serial == ndna:
       sc_index = bead_tot_dna + 2*residue_id + 1
   else:
       count = bead_tot_dna
       for i in range (ndna, serial):
          count = count + info_dna_pro[name_list[i]]
       sc_index = count + 2*residue_id + 1
   return sc_index

#Note: Get bead index of SOP-SC backbone bead

def get_pro_bb_index(serial, residue_id):
   if serial == ndna:
       bb_index = bead_tot_dna + 2*residue_id
   else:
       count = bead_tot_dna
       for i in range (ndna, serial):
          count = count + info_dna_pro[name_list[i]]
       bb_index = count + 2*residue_id
   return bb_index

#Note: Calculates sum total of beads before the current bead 

def get_sum_pre_pro_bead_number(serial):
   if serial == ndna:
      bead_no = bead_tot_dna
   else:
      count = bead_tot_dna
      for i in range (ndna, serial):
          count = count + info_dna_pro[name_list[i]]
      bead_no = count
   return bead_no

#Note: Get radius of SOP-SC sidechain bead

def get_radii_pro(serial, residue_id):
   if serial == ndna:
       sc_index = bead_tot_dna + 2*residue_id + 1
       key = residue[sc_index]
       radius = radii_pro[key]
   else:
       count = bead_tot_dna
       for i in range (ndna, serial):
          count = count + info_dna_pro[name_list[i]]
       sc_index = count + 2*residue_id + 1
       key = residue[sc_index]
       radius = radii_pro[key]
   return radius

#Note: Checks if a bead is odd or even

def even(ind):
    if (ind % 2) == 0:
       return True
    else:
       return False

#Note: Bond distance between two SOP-SC beads

def get_dist_along_chain(ind1, ind2):

    #NOTE: This function is applicable in case of dsDNA, modify if the DNA is ss

    bead_ind1 = ind1 - bead_tot_dna
    bead_ind2 = ind2 - bead_tot_dna

    if even(bead_ind1) == False and even(bead_ind2) == False:
       dist = abs((bead_ind1 - bead_ind2)/2) + 2
       return int(dist)
    if even(bead_ind1) == True and even(bead_ind2) == True:
       dist = abs((bead_ind1 - bead_ind2)/2)
       return int(dist)
    if even(bead_ind1) == True and even(bead_ind2) == False:
       dist = abs((bead_ind1 - (bead_ind2-1))/2) + 1
       return int(dist)
    if even(bead_ind1) == False and even(bead_ind2) == True:
       dist = abs(((bead_ind1-1) - bead_ind2)/2) + 1
       return int(dist)
#Note: Computes interaction strngth from Betancourt Thirumalai statistical potential

def get_eps_from_BT_pot(resid1, resid2):
    index1 = Pro_Seq_BT.index(resid1)
    index2 = Pro_Seq_BT.index(resid2)

    if index2 > index1:
        eps = BT_St_Pot[index1][index2]
    else:
        eps = BT_St_Pot[index2][index1]

    return eps

