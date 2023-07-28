import numpy as np
import math
from user_input import *
from collections import Counter


# Constants:

k_B = 0.001987204259    #Boltzmann Constant in kcal/(mol-K)


# Electrostatics

diele_pro = 20.0         	#Dielectric constant for electrostatic interactions for Proteins
diele_dna = 78.0         	#Dielectric constant for electrostatic interactions for DNA
diele_dna_pro = 78.0     	#Dielectric constant for electrostatic interactions between DNA and Proteins
charge_per_uni_len = 4.4        #Charge per unit length of DNA (dsDNA) in Angstrom
bjerrum_length = (2.56404*10000000)/(diele_dna*temp*153.616394935) #Bjerrum length in Angstrom
norm_charge_p = -(charge_per_uni_len/bjerrum_length)      #Renormalized charge on phosphate group following Manning Condensation (see main paper)
coulomb_prefix_val = 332.0637090  # e*e*NA/4.0*pi*eps_0 (e = elementary charge, NA = Avogadro's constant, eps_0 = vacuum permittivity)
number_density = 2.0*salt_conc*6.022*0.0000001
kappa_dna_val = np.sqrt((4 * np.pi * (coulomb_prefix_val / diele_dna)* number_density)/(k_B * temp))
kappa_pro_val = np.sqrt((4 * np.pi * (coulomb_prefix_val / diele_pro)* number_density)/(k_B * temp))
kappa_dna_pro_val = np.sqrt((4 * np.pi * (coulomb_prefix_val / diele_dna_pro)* number_density)/(k_B * temp))

# Parameters for SOP-SC model

bond_sep_nat = 3                #Minimum bond separation for backbone-backbone (BB-BB) and backbone-sidechain (BB-SC) native contact
bond_sep_nat_sc = 4             #Minimum bond separation for SC-SC native contact
radii_scale = 0.8                  #Scaling factor applicable to angle potentials between BB and SC beads 
pro_bb_radius = 1.9             #Radius of BB beads in Angstrom
pro_nat_cut = 8.00              #Protein beads within this cut-off distance (in angstrom) in PDB will have native interactions 
pro_dna_nat_cut = 11.00         #Cut-off distance in Angstrom between DNA and protein beads for native interaction
pro_eps_bb = 0.5                #Native interaction strength between BB beads (in kcal/mol) 
pro_eps_bs = 0.5                #Native interaction strength between BB and SC beads (in kcal/mol)
pro_eps_ss = 0.3                #Native interaction strength between SC beads (in kcal/mol)
kfene = 20.0                    #Force constant for FENE potential in kcal/mol/angstrom**2
eps_angle_pro = 1.0             #Interaction strength for angle potential for Proteins in kcal/mol
R_max_fene = 2.0                #Tolerance in Fene potential in angstroms
scale_nat_dna_pro = 0.15        #Scale factor for DNA-Protein native interaction


# Parameters for TIS-DNA model

# Base Stacking

kbond_st = 1.45               #Force constant for bonded interaction in stacking in (Angstrom*Angstrom)**-1
kphi1_st = 3.0                #Force constant for angular interaction in stacking in (radians*radians)**-1
kphi2_st = 3.0                #Force constant for torsional interaction in stacking in (radians*radians)**-1

# Hydrogen Bonding

en_hbond = -1.92              #Hydrogen bonding energies in kcal/mol 
pi_cut = 3.124139361          #Angle values close to pi to avoid instability
kbond_hb = 5.00               #Force constant for bonded interaction in hydrogen bonding in (Angstrom*Angstrom)**-1        
kang_hb = 1.5                 #Force constant for angular interaction in hydrogen bonding in (radians*radians)**-1
kdih_hb = 0.15                #Force constant for dihedral interaction in hydrogen bonding in (radians*radians)**-1
