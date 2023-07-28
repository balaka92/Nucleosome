from global_param import *
import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~Matrix denoting effective energy scales between P of DNA and amino acid BB~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

P_BB = {"ALA" : {"DG" : 0.41, "DA" : 0.27, "DC" : 0.31, "DT" : 0.40}, "GLY" : {"DG" : 0.62, "DA" : 0.58, "DC" : 0.52, "DT" : 0.54}, "LEU" : {"DG" : 0.17, "DA" : 0.08, "DC" : 0.15, "DT" : 0.19}, "ILE" : {"DG" : 0.38, "DA" : 0.34, "DC" : 0.33, "DT" : 0.29}, "TYR" : {"DG" : 0.33, "DA" : 0.24, "DC" : 0.32, "DT" : 0.30}, "TRP" : {"DG" : 0.27, "DA" : 0.33, "DC" : 0.35, "DT" : 0.35}, "LYS" : {"DG" : 0.69, "DA" : 0.62, "DC" : 0.66, "DT" : 0.66}, "PRO" : {"DG" : 0.37, "DA" : 0.36, "DC" : 0.35, "DT" : 0.39}, "GLU" : {"DG" : 0.12, "DA" : 0.06, "DC" : 0.13, "DT" : 0.04}, "HIS" : {"DG" : 0.72, "DA" : 0.58, "DC" : 0.43, "DT" : 0.56}, "ASP" : {"DG" : 0.38, "DA" : 0.20, "DC" : 0.31, "DT" : 0.01}, "ARG" : {"DG" : 0.82, "DA" : 0.68, "DC" : 0.65, "DT" : 0.70}, "GLN" : {"DG" : 0.49, "DA" : 0.41, "DC" : 0.41, "DT" : 0.43}, "ASN" : {"DG" : 0.60, "DA" : 0.51, "DC" : 0.50, "DT" : 0.55}, "CYS" : {"DG" : 0.54, "DA" : 0.22, "DC" : 0.25, "DT" : 0.27}, "SER" : {"DG" : 0.70, "DA" : 0.56, "DC" : 0.59, "DT" : 0.60}, "MET" : {"DG" : 0.36, "DA" : 0.33, "DC" : 0.34, "DT" : 0.30}, "VAL" : {"DG" : 0.43, "DA" : 0.32, "DC" : 0.25, "DT" : 0.32}, "THR" : {"DG" : 0.67, "DA" : 0.57, "DC" : 0.57, "DT" : 0.64}, "PHE" : {"DG" : 0.44, "DA" : 0.25, "DC" : 0.26, "DT" : 0.33}}



#~~~~~~~~~~~~~~~~~~~~~~~~~Matrix denoting effective energy scales between S of DNA and amino acid BB~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

S_BB = {"ALA" : {"DG" : 0.43, "DA" : 0.30, "DC" : 0.38, "DT" : 0.34}, "GLY" : {"DG" : 0.71, "DA" : 0.69, "DC" : 0.66, "DT" : 0.62}, "LEU" : {"DG" : 0.07, "DA" : 0.03, "DC" : 0.06, "DT" : 0.06}, "ILE" : {"DG" : 0.33, "DA" : 0.31, "DC" : 0.30, "DT" : 0.22}, "TYR" : {"DG" : 0.33, "DA" : 0.30, "DC" : 0.35, "DT" : 0.24}, "TRP" : {"DG" : 0.20, "DA" : 0.44, "DC" : 0.20, "DT" : 0.32}, "LYS" : {"DG" : 0.74, "DA" : 0.66, "DC" : 0.71, "DT" : 0.66}, "PRO" : {"DG" : 0.42, "DA" : 0.43, "DC" : 0.42, "DT" : 0.41}, "GLU" : {"DG" : 0.11, "DA" : 0.02, "DC" : 0.07, "DT" : 0.00}, "HIS" : {"DG" : 0.65, "DA" : 0.66, "DC" : 0.51, "DT" : 0.57}, "ASP" : {"DG" : 0.37, "DA" : 0.24, "DC" : 0.37, "DT" : 0.02}, "ARG" : {"DG" : 0.83, "DA" : 0.75, "DC" : 0.73, "DT" : 0.75}, "GLN" : {"DG" : 0.51, "DA" : 0.48, "DC" : 0.48, "DT" : 0.40}, "ASN" : {"DG" : 0.69, "DA" : 0.61, "DC" : 0.60, "DT" : 0.63}, "CYS" : {"DG" : 0.38, "DA" : 0.16, "DC" : 0.24, "DT" : 0.28}, "SER" : {"DG" : 0.85, "DA" : 0.63, "DC" : 0.74, "DT" : 0.66}, "MET" : {"DG" : 0.25, "DA" : 0.31, "DC" : 0.27, "DT" : 0.24}, "VAL" : {"DG" : 0.42, "DA" : 0.32, "DC" : 0.25, "DT" : 0.31}, "THR" : {"DG" : 0.71, "DA" : 0.66, "DC" : 0.62, "DT" : 0.70}, "PHE" : {"DG" : 0.40, "DA" : 0.20, "DC" : 0.33, "DT" : 0.30}}


#~~~~~~~~~~~~~~~~~~~~~~~~~Matrix denoting effective energy scales between S of DNA and amino acid SC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

S_SC = {"ALA" : {"DG" : 0.56, "DA" : 0.44, "DC" : 0.47, "DT" : 0.45}, "GLY" : {"DG" : 0.80, "DA" : 0.79, "DC" : 0.73, "DT" : 0.73}, "LEU" : {"DG" : 0.06, "DA" : 0.07, "DC" : 0.01, "DT" : 0.11}, "ILE" : {"DG" : 0.31, "DA" : 0.33, "DC" : 0.25, "DT" : 0.24}, "TYR" : {"DG" : 0.62, "DA" : 0.58, "DC" : 0.61, "DT" : 0.57}, "TRP" : {"DG" : 0.38, "DA" : 0.45, "DC" : 0.45, "DT" : 0.41}, "LYS" : {"DG" : 0.97, "DA" : 0.90, "DC" : 0.89, "DT" : 0.89}, "PRO" : {"DG" : 0.49, "DA" : 0.50, "DC" : 0.49, "DT" : 0.50}, "GLU" : {"DG" : 0.32, "DA" : 0.26, "DC" : 0.31, "DT" : 0.21}, "HIS" : {"DG" : 0.83, "DA" : 0.86, "DC" : 0.75, "DT" : 0.78}, "ASP" : {"DG" : 0.54, "DA" : 0.39, "DC" : 0.56, "DT" : 0.18}, "ARG" : {"DG" : 1.14, "DA" : 0.99, "DC" : 0.94, "DT" : 1.02}, "GLN" : {"DG" : 0.68, "DA" : 0.69, "DC" : 0.67, "DT" : 0.55}, "ASN" : {"DG" : 0.81, "DA" : 0.75, "DC" : 0.76, "DT" : 0.79}, "CYS" : {"DG" : 0.56, "DA" : 0.28, "DC" : 0.28, "DT" : 0.42}, "SER" : {"DG" : 0.94, "DA" : 0.75, "DC" : 0.86, "DT" : 0.78}, "MET" : {"DG" : 0.35, "DA" : 0.38, "DC" : 0.45, "DT" : 0.30}, "VAL" : {"DG" : 0.42, "DA" : 0.39, "DC" : 0.33, "DT" : 0.36}, "THR" : {"DG" : 0.82, "DA" : 0.78, "DC" : 0.76, "DT" : 0.85}, "PHE" : {"DG" : 0.47, "DA" : 0.36, "DC" : 0.41, "DT" : 0.37}}


#~~~~~~~~~~~~~~~~~~~~~~~~~Matrix denoting effective energy scales between B of DNA and amino acid BB~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

B_BB = {"ALA" : {"DG" : 0.69, "DA" : 0.68, "DC" : 0.67, "DT" : 0.81}, "GLY" : {"DG" : 1.11, "DA" : 1.04, "DC" : 1.11, "DT" : 1.04}, "LEU" : {"DG" : 0.00, "DA" : 0.10, "DC" : 0.14, "DT" : 0.35}, "ILE" : {"DG" : 0.37, "DA" : 0.44, "DC" : 0.47, "DT" : 0.52}, "TYR" : {"DG" : 0.60, "DA" : 0.51, "DC" : 0.72, "DT" : 0.58}, "TRP" : {"DG" : 0.33, "DA" : 0.60, "DC" : 0.58, "DT" : 0.62}, "LYS" : {"DG" : 0.96, "DA" : 0.83, "DC" : 1.01, "DT" : 0.85}, "PRO" : {"DG" : 0.58, "DA" : 0.62, "DC" : 0.61, "DT" : 0.67}, "GLU" : {"DG" : 0.51, "DA" : 0.22, "DC" : 0.48, "DT" : 0.30}, "HIS" : {"DG" : 1.11, "DA" : 0.96, "DC" : 1.00, "DT" : 1.04}, "ASP" : {"DG" : 0.86, "DA" : 0.45, "DC" : 0.86, "DT" : 0.37}, "ARG" : {"DG" : 1.18, "DA" : 1.02, "DC" : 1.16, "DT" : 1.13}, "GLN" : {"DG" : 0.86, "DA" : 0.94, "DC" : 0.92, "DT" : 0.98}, "ASN" : {"DG" : 1.11, "DA" : 1.09, "DC" : 1.07, "DT" : 1.19}, "CYS" : {"DG" : 0.62, "DA" : 0.60, "DC" : 0.70, "DT" : 0.58}, "SER" : {"DG" : 1.16, "DA" : 1.02, "DC" : 1.18, "DT" : 1.08}, "MET" : {"DG" : 0.38, "DA" : 0.59, "DC" : 0.66, "DT" : 0.68}, "VAL" : {"DG" : 0.62, "DA" : 0.58, "DC" : 0.56, "DT" : 0.65}, "THR" : {"DG" : 0.93, "DA" : 0.95, "DC" : 0.97, "DT" : 1.03}, "PHE" : {"DG" : 0.61, "DA" : 0.43, "DC" : 0.62, "DT" : 0.65}}


#~~~~~~~~~~~~~~~~~~~~~~~~~Matrix denoting effective energy scales between B of DNA and amino acid SC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

B_SC = {"ALA" : {"DG" : 0.75, "DA" : 0.71, "DC" : 0.70, "DT" : 0.78}, "GLY" : {"DG" : 1.07, "DA" : 1.02, "DC" : 1.04, "DT" : 1.01}, "LEU" : {"DG" : 0.00, "DA" : 0.03, "DC" : 0.01, "DT" : 0.24}, "ILE" : {"DG" : 0.32, "DA" : 0.38, "DC" : 0.43, "DT" : 0.40}, "TYR" : {"DG" : 0.86, "DA" : 0.74, "DC" : 0.90, "DT" : 0.85}, "TRP" : {"DG" : 0.21, "DA" : 0.20, "DC" : 0.46, "DT" : 0.39}, "LYS" : {"DG" : 1.13, "DA" : 1.09, "DC" : 1.22, "DT" : 1.12}, "PRO" : {"DG" : 0.63, "DA" : 0.64, "DC" : 0.58, "DT" : 0.69}, "GLU" : {"DG" : 0.69, "DA" : 0.42, "DC" : 0.66, "DT" : 0.43}, "HIS" : {"DG" : 1.21, "DA" : 1.13, "DC" : 1.18, "DT" : 1.17}, "ASP" : {"DG" : 1.00, "DA" : 0.59, "DC" : 0.98, "DT" : 0.53}, "ARG" : {"DG" : 1.52, "DA" : 1.38, "DC" : 1.49, "DT" : 1.41}, "GLN" : {"DG" : 1.09, "DA" : 1.06, "DC" : 1.03, "DT" : 1.08}, "ASN" : {"DG" : 1.22, "DA" : 1.20, "DC" : 1.21, "DT" : 1.30}, "CYS" : {"DG" : 0.56, "DA" : 0.50, "DC" : 0.68, "DT" : 0.44}, "SER" : {"DG" : 1.27, "DA" : 1.11, "DC" : 1.24, "DT" : 1.14}, "MET" : {"DG" : 0.51, "DA" : 0.62, "DC" : 0.63, "DT" : 0.72}, "VAL" : {"DG" : 0.63, "DA" : 0.61, "DC" : 0.55, "DT" : 0.64}, "THR" : {"DG" : 1.09, "DA" : 1.07, "DC" : 1.08, "DT" : 1.15}, "PHE" : {"DG" : 0.64, "DA" : 0.50, "DC" : 0.54, "DT" : 0.62}}


#~~~~~~~~Format is mass, radius, charge, eps~~~~~~~~~~~~~~~~#

dna_bead = {"DA": (130.0335,2.8,0.0,0.2), "DT": (120.0014,2.7,0.0,0.2), "DC": (106.0141,2.7,0.0,0.2), "DG": (146.0275,3.0,0.0,0.2), "DS2" : (75.994, 2.9, 0.0, 0.2), "DP2": (94.9498, 2.0,-1.0,0.2)}

#~~~~~~~Radii of SOP-SC beads~~~~~~~~~#

radii_pro = {"GLY": 0.5, "ALA": 2.52, "VAL": 2.93, "LEU": 3.09, "ILE": 3.09, "MET": 3.09, "PHE" : 3.18, "PRO" : 2.78, "SER": 2.59, "THR" : 2.81, "ASN" : 2.84, "GLN" : 3.01, "TYR": 3.23, "TRP" : 3.39, "ASP": 2.79, "GLU":2.96, "HIS": 3.04, "LYS": 3.18, "ARG": 3.28, "CYS" : 2.74}


#~~~~~~~~~Force constant for different bonds in TIS-DNA~~~~~~~~~~~~#
# values in kcal/mol/angstrom**2

bond_type = {"DS2DP2" : 62.59, "DP2DS2" : 17.63, "DS2DA" : 44.31, "DS2DT" : 46.56, "DS2DC" : 43.25, "DS2DG": 48.98}


#~~~~~~~~~Force constant for different angles in TIS-DNA~~~~~~~~~~~~#
# values in kcal/mol/radian**2


angle_type = {"DP2DS2DP2" : 25.67, "DS2DP2DS2" : 67.50, "DP2DS2DA" : 29.53, "DP2DS2DT" : 39.56, "DP2DS2DG" : 26.28, "DP2DS2DC" : 35.02, "DADS2DP2" : 67.32, "DTDS2DP2" : 93.99, "DGDS2DP2" : 62.94, "DCDS2DP2" : 77.78}


#~~~~~~~~~~~Stacking parameter~~~~~~~~~~~~#

stack_type = {"DADA": (5.69, 0.94, 322.0), "DADT" : (4.95, 0.87, 293.0), "DTDA" : (5.02, 0.65, 293.0), "DADC" : (4.98, 0.92, 293.0), "DCDA" : (5.03, 0.78, 293.0), "DADG" : (5.43, 1.06, 333.6), "DGDA" : (5.42, 0.98, 333.6), "DCDT" : (4.13, 0.71, 288.3), "DTDC" : (4.14, 0.94, 288.3), "DCDC" : (4.12, 0.92, 288.3), "DCDG" : (5.21, 2.31, 331.9), "DGDC" : (5.13, 2.41, 331.9), "DGDT" : (5.28, 1.69, 332.6), "DTDG" : (5.43, 1.72, 332.6), "DGDG" : (5.66, -0.29, 353.9), "DTDT" : (4.17, 0.89, 288.3)}


#~~~~~~~~~~~~Strength of Hydrogen bonding~~~~~~~~~~~~~~#

hbond_type = {"DADT" : 2.0*en_hbond, "DTDA" : 2.0*en_hbond, "DGDC" : 3.0*en_hbond, "DCDG" : 3.0*en_hbond}

#~~~~~~~~~~~~Betancourt Thirumalai Statistical Potential~~~~~~~~~~~~~~~#

BT_filename = "BT_Col.csv" #File containing upper triangular part of the BT matrix
BT_St_Pot = np.genfromtxt(BT_filename, delimiter=',')    # Reading Betancourt Thirumalai Statistical Potential #
Pro_Seq_BT = ["CYS","PHE","LEU","TRP","VAL","ILE","MET","HIS","TYR","ALA","GLY","PRO","ASN","THR","SER","ARG","GLN","ASP","LYS","GLU"]
