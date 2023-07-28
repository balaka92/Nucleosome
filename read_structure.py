import numpy as np
from global_param import *
from collections import Counter

atom_index = list()
residue = list()
atom_type = list()
chain_ID = list()
residue_no = list()
coor_x = list()
coor_y = list()
coor_z = list()

def read_pdb(filename):
    for line in open(filename):
       id=str(line[0:6]).strip()
       if id == 'ATOM':
          atom_index.append(int(line[6:11]))
          atom_type.append(str(line[12:16]).strip())
          residue.append(str(line[17:20]).strip())
          chain_ID.append(str(line[21:22]))
          residue_no.append(int(line[23:26]))
          coor_x.append(float(line[30:38]))
          coor_y.append(float(line[39:46]))
          coor_z.append(float(line[47:54]))
    return atom_index, atom_type, residue, chain_ID, residue_no, coor_x, coor_y, coor_z
atom_index, atom_type, residue, chain_ID, residue_no, coor_x, coor_y, coor_z = read_pdb(pdb_filename)

dna_list = []
pro_list = []
for i in range(0,len(atom_index)):
    atm_type = atom_type[i]
    if atm_type == "DP2" or atm_type == "DS2" or atm_type == "DB2":
        dna_list.append(chain_ID[i])
    else:
        pro_list.append(chain_ID[i])
ndna = len(Counter(dna_list).keys())
npro = len(Counter(pro_list).keys())
nchain=ndna+npro
info_dna_pro = dict(Counter(chain_ID))
name_list = list(info_dna_pro.keys())
dict_nres = {}
dict_nbead = {}
count_dna=0
count_pro=0
for i in range(len(info_dna_pro)):
    key = list(info_dna_pro.keys())[i]
    name="nres_tot_chain_" + str(key)
    name2="bead_tot_chain_" + str(key)
    if i < ndna:
        val=int((info_dna_pro[key] + 1)/3)
        count_dna=count_dna + int(info_dna_pro[key])
    else:
        val=int((info_dna_pro[key])/2)
        count_pro=count_pro + int(info_dna_pro[key])
    dict_nres[name] = val
    dict_nbead[name2] = int(info_dna_pro[key])

bead_tot_dna=count_dna
bead_tot_pro=count_pro
bead_tot=bead_tot_dna+bead_tot_pro

