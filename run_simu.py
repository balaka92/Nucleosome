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
from create_topology import *

#~~~~~~~~~Custom Energy Reporter~~~~~~~~~~#

class EnergyCompReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, None)

    def report(self, simulation, state):

        openmm_state = simulation.context.getState(getEnergy=True)
        potential_energy = openmm_state.getPotentialEnergy()
        bond_dna = simulation.context.getState(getEnergy=True,groups={0}).getPotentialEnergy()
        bond_pro = simulation.context.getState(getEnergy=True,groups={1}).getPotentialEnergy()
        ang_dna = simulation.context.getState(getEnergy=True,groups={2}).getPotentialEnergy()
        ang_pro = simulation.context.getState(getEnergy=True,groups={3}).getPotentialEnergy()
        st_dna = simulation.context.getState(getEnergy=True,groups={4}).getPotentialEnergy()
        hb_dna = simulation.context.getState(getEnergy=True,groups={5}).getPotentialEnergy()
        nat_pro = simulation.context.getState(getEnergy=True,groups={6}).getPotentialEnergy()
        nat_pro_dna = simulation.context.getState(getEnergy=True,groups={7}).getPotentialEnergy()
        ev_dna_dna = simulation.context.getState(getEnergy=True,groups={8}).getPotentialEnergy()
        ev_pro_pro = simulation.context.getState(getEnergy=True,groups={9}).getPotentialEnergy()
        ev_dna_pro_P_SC = simulation.context.getState(getEnergy=True,groups={10}).getPotentialEnergy()
        ev_dna_pro_P_BB = simulation.context.getState(getEnergy=True,groups={11}).getPotentialEnergy()
        ev_dna_pro_rest = simulation.context.getState(getEnergy=True,groups={12}).getPotentialEnergy()
        el_dna_dna = simulation.context.getState(getEnergy=True,groups={13}).getPotentialEnergy()
        el_pro_pro = simulation.context.getState(getEnergy=True,groups={14}).getPotentialEnergy()
        el_dna_pro = simulation.context.getState(getEnergy=True,groups={15}).getPotentialEnergy()


        self._out.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s \t%s\n" % (simulation.currentStep,potential_energy/potential_energy.unit, bond_dna/bond_dna.unit, bond_pro/bond_pro.unit, ang_dna/ang_dna.unit, ang_pro/ang_pro.unit, st_dna/st_dna.unit, hb_dna/hb_dna.unit, nat_pro/nat_pro.unit, nat_pro_dna/nat_pro_dna.unit, ev_dna_dna/ev_dna_dna.unit, ev_pro_pro/ev_pro_pro.unit, ev_dna_pro_P_SC/ev_dna_pro_P_SC.unit, ev_dna_pro_P_BB/ev_dna_pro_P_BB.unit, ev_dna_pro_rest/ev_dna_pro_rest.unit, el_dna_dna/el_dna_dna.unit, el_pro_pro/el_pro_pro.unit, el_dna_pro/el_dna_pro.unit))




#////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Simulation Block~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#



time_st = datetime.datetime.utcnow()

pdb = app.PDBFile(pdb_filename)
positions = pdb.positions
system = mm.System()

for i in range (0, bead_tot):
    mass_i = 50.0
    system.addParticle(mass_i*u.amu)

################################################Bonds######################################################

###############DNA BONDS###################

F_Bond_DNA = mm.HarmonicBondForce()

F_Bond_DNA.setUsesPeriodicBoundaryConditions(False)

for i in range (0, len(bond_dna)):
    particle1 = int(bond_dna[i][0])
    particle2 = int(bond_dna[i][1])
    bond_k = 2.0*float(bond_dna[i][2])
    bond_r0 = float(bond_dna[i][3])

    F_Bond_DNA.addBond(particle1, particle2, bond_r0*u.angstroms, bond_k*u.kilocalories_per_mole/u.angstroms**2)



###########PRO BONDS###########

F_Bond_PRO = mm.CustomCompoundBondForce(2,"-0.5*kfene*R_max_fene*R_max_fene*log(1.0 - (((distance(p1,p2) - rbond)*(distance(p1,p2) - rbond))/(R_max_fene*R_max_fene)));")

F_Bond_PRO.addGlobalParameter("kfene", kfene*u.kilocalories_per_mole/u.angstroms**2)
F_Bond_PRO.addGlobalParameter("R_max_fene", R_max_fene*u.angstroms)
F_Bond_PRO.addPerBondParameter("rbond")

F_Bond_PRO.setUsesPeriodicBoundaryConditions(False)


for i in range (0, len(bond_pro)):
    p1 = int(bond_pro[i][0])
    p2 = int(bond_pro[i][1])
    rbond = float(bond_pro[i][2])
    group_bond = [p1,p2]
    F_Bond_PRO.addBond(group_bond, [rbond*u.angstroms])



################################################Angles######################################################


##########DNA ANGLES############


F_Angle_DNA = mm.HarmonicAngleForce()

F_Angle_DNA.setUsesPeriodicBoundaryConditions(False)


for i in range (0, len(angle_dna)):
    particle1 = int(angle_dna[i][0])
    particle2 = int(angle_dna[i][1])
    particle3 = int(angle_dna[i][2])
    ang_k = 2.0*float(angle_dna[i][3])
    ang_theta = float(angle_dna[i][4])

    F_Angle_DNA.addAngle(particle1, particle2, particle3, ang_theta*u.radians, ang_k*u.kilocalories_per_mole/u.radians**2)




#########PRO ANGLES##############

F_Angle_PRO = mm.CustomCompoundBondForce(2,"eps_angle_pro*(sigma/(distance(p1,p2)))^6;")

F_Angle_PRO.addGlobalParameter("eps_angle_pro", eps_angle_pro*u.kilocalories_per_mole)
F_Angle_PRO.addPerBondParameter("sigma")

F_Angle_PRO.setUsesPeriodicBoundaryConditions(False)


for i in range (0, len(angle_pro)):
    p1 = int(angle_pro[i][0])
    p2 = int(angle_pro[i][1])
    sigma = float(angle_pro[i][2])
    group_angle = [p1,p2]
    F_Angle_PRO.addBond(group_angle, [sigma*u.angstroms])



################~~~~~~~~~~~Stacking within DNA~~~~~~~~~~~####################


Stackingforce = mm.CustomCompoundBondForce(7, "U0/(1.0 + kbond_st*(distance(p6, p7) - r0)^2  + kphi1_st*(dihedral(p2,p1,p3,p4) - phi10 + pi*(((pi - dihedral(p2,p1,p3,p4) + phi10)/abs(pi - dihedral(p2,p1,p3,p4) + phi10)) -  ((pi + dihedral(p2,p1,p3,p4) - phi10)/abs(pi + dihedral(p2,p1,p3,p4) - phi10)))  )^2  +  kphi2_st*(dihedral(p5,p4,p3,p1) - phi20 + pi*(((pi - dihedral(p5,p4,p3,p1) + phi20)/abs(pi - dihedral(p5,p4,p3,p1) + phi20)) -  ((pi + dihedral(p5,p4,p3,p1) - phi20)/abs(pi + dihedral(p5,p4,p3,p1) - phi20)))  )^2 )");
Stackingforce.addPerBondParameter("U0");
Stackingforce.addPerBondParameter("r0");
Stackingforce.addPerBondParameter("phi10");
Stackingforce.addPerBondParameter("phi20");
Stackingforce.addGlobalParameter('kbond_st', kbond_st/u.angstroms**2)
Stackingforce.addGlobalParameter('kphi1_st', kphi1_st/u.radians**2)
Stackingforce.addGlobalParameter('kphi2_st', kphi2_st/u.radians**2)
Stackingforce.addGlobalParameter('pi', np.pi)
Stackingforce.setUsesPeriodicBoundaryConditions(False)


for i in range (0, len(stacked)):
    p1 = int(stacked[i][0])
    p2 = int(stacked[i][1])
    p3 = int(stacked[i][2])
    p4 = int(stacked[i][3])
    p5 = int(stacked[i][4])
    p6 = int(stacked[i][5])
    p7 = int(stacked[i][6])
    U0 = float(stacked[i][7])
    r0 = float(stacked[i][8])
    phi10 = float(stacked[i][9])
    phi20 = float(stacked[i][10])
    group_st = [p1,p2,p3,p4,p5,p6,p7]

    Stackingforce.addBond(group_st, [U0*u.kilocalories_per_mole,r0*u.angstroms, phi10*u.radians, phi20*u.radians])




###############~~~~~~~~~~~~Hydrogen Bonding within DNA~~~~~~~~~~~#################
#NOTE: Do not compute force if angles approach 180 degrees


Hbondingforce = mm.CustomCompoundBondForce(6,"scale1*UHyd/(1.0 + kbond_hb*DIST*DIST + kang_hb*ANG1*ANG1 + kang_hb*ANG2*ANG2 + kdih_hb*DIHED1*DIHED1 + kdih_hb*DIHED2*DIHED2 + kdih_hb*DIHED3*DIHED3); scale1 = step(pi_cut - angle(p2, p1, p4))*step(pi_cut - angle(p5, p4, p1)); DIST = (distance(p1,p4) - dist_0); ANG1 = (angle(p2,p1,p4) - angle_01); ANG2 = (angle(p5,p4,p1) - angle_02); DIHED1 = dihedral(p2,p1,p4,p5) - dihedral_01 + pi*(((pi - dihedral(p2,p1,p4,p5) + dihedral_01)/abs(pi - dihedral(p2,p1,p4,p5) + dihedral_01)) - ((pi + dihedral(p2,p1,p4,p5) - dihedral_01)/abs(pi + dihedral(p2,p1,p4,p5) - dihedral_01))) ;  DIHED2 = (dihedral(p6,p5,p4,p1) - dihedral_02 + pi*(((pi - dihedral(p6,p5,p4,p1) + dihedral_02)/abs(pi - dihedral(p6,p5,p4,p1) + dihedral_02)) - ((pi + dihedral(p6,p5,p4,p1) - dihedral_02)/abs(pi + dihedral(p6,p5,p4,p1) - dihedral_02)))); DIHED3 = (dihedral(p3,p2,p1,p4) - dihedral_03 + pi*(((pi - dihedral(p3,p2,p1,p4) + dihedral_03)/abs(pi - dihedral(p3,p2,p1,p4) + dihedral_03)) - ((pi + dihedral(p3,p2,p1,p4) - dihedral_03)/abs(pi + dihedral(p3,p2,p1,p4) - dihedral_03))));")

Hbondingforce.addPerBondParameter("UHyd")
Hbondingforce.addPerBondParameter("dist_0")
Hbondingforce.addPerBondParameter("angle_01")
Hbondingforce.addPerBondParameter("angle_02")
Hbondingforce.addPerBondParameter("dihedral_01")
Hbondingforce.addPerBondParameter("dihedral_02")
Hbondingforce.addPerBondParameter("dihedral_03")
Hbondingforce.addGlobalParameter('pi_cut',pi_cut*u.radians)
Hbondingforce.addGlobalParameter('kbond_hb', kbond_hb/u.angstroms**2)
Hbondingforce.addGlobalParameter('kang_hb', kang_hb/u.radians**2)
Hbondingforce.addGlobalParameter('kdih_hb', kdih_hb/u.radians**2)
Hbondingforce.addGlobalParameter('pi', np.pi)
Hbondingforce.setUsesPeriodicBoundaryConditions(False)



for i in range (0, len(hbonded)):
    p1 = int(hbonded[i][0])
    p2 = int(hbonded[i][1])
    p3 = int(hbonded[i][2])
    p4 = int(hbonded[i][3])
    p5 = int(hbonded[i][4])
    p6 = int(hbonded[i][5])
    dist_0 = float(hbonded[i][6])
    angle_01 = float(hbonded[i][7])
    angle_02 = float(hbonded[i][8])
    dihedral_01 = float(hbonded[i][9])
    dihedral_02 = float(hbonded[i][10])
    dihedral_03 = float(hbonded[i][11])
    UHyd = float(hbonded[i][12])
    group_hb = [p1,p2,p3,p4,p5,p6]

    Hbondingforce.addBond(group_hb, [UHyd*u.kilocalories_per_mole, dist_0*u.angstroms, angle_01*u.radians, angle_02*u.radians, dihedral_01*u.radians, dihedral_02*u.radians, dihedral_03*u.radians])
#system.addForce(Hbondingforce)



############~~~~~~~~~~Native Interactions Within PRO~~~~~~~~~~################

F_Native_PRO = mm.CustomCompoundBondForce(2,"eps_nat*(((rcry/distance(p1,p2))^12) - 2.0*((rcry/(distance(p1,p2)))^6))")
F_Native_PRO.addPerBondParameter("rcry")
F_Native_PRO.addPerBondParameter("eps_nat")
F_Native_PRO.setUsesPeriodicBoundaryConditions(False)



for i in range (0, len(native_pro)):
    p1 = int(native_pro[i][0])
    p2 = int(native_pro[i][1])
    rcry = float(native_pro[i][2])
    eps_nat = float(native_pro[i][3])
    group_n_p = [p1,p2]

    F_Native_PRO.addBond(group_n_p, [rcry*u.angstroms,eps_nat*u.kilocalories_per_mole])


############~~~~~~~~~~~Native Interactions between PRO and DNA~~~~~~~~~~~~############

F_Native_DNA_PRO = mm.CustomCompoundBondForce(2, "eps_dna_pro_nat*(((rcry_DP/(distance(p1,p2)))^12) - (2.0*(rcry_DP/(distance(p1,p2)))^6))")

F_Native_DNA_PRO.addPerBondParameter("eps_dna_pro_nat")
F_Native_DNA_PRO.addPerBondParameter("rcry_DP")
F_Native_DNA_PRO.setUsesPeriodicBoundaryConditions(False)

for i in range (0, len(dna_pro_native)):
    p1 = int(dna_pro_native[i][0])
    p2 = int(dna_pro_native[i][1])
    rcry_DP = float(dna_pro_native[i][2])
    eps_dna_pro_nat = float(dna_pro_native[i][3])
    group_n_d_p = [p1,p2]
    F_Native_DNA_PRO.addBond(group_n_d_p, [eps_dna_pro_nat*u.kilocalories_per_mole, rcry_DP*u.angstroms])




#~~~~~~~~~~~~~~~~~~~Excluded Volume interactions~~~~~~~~~~~~~~~~~~~#


# dna-dna

F_EV_DNA_DNA = mm.CustomNonbondedForce( 'step(sig-r)*eps*((sig/r)^12 - 2.0*(sig/r)^6 + 1.0); eps=sqrt(eps1*eps2);')
F_EV_DNA_DNA.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

F_EV_DNA_DNA.addPerParticleParameter('eps')
F_EV_DNA_DNA.addGlobalParameter("sig", 3.2*u.angstroms)

for i in range (0, len(dna_dna_ex_vol)):
    eps = float(dna_dna_ex_vol[i][0])

    F_EV_DNA_DNA.addParticle([eps*u.kilocalories_per_mole])
for j in range (0, len(exclusion_list)):
    particle1 = int(exclusion_list[j][0])
    particle2 = int(exclusion_list[j][1])

    F_EV_DNA_DNA.addExclusion(particle1, particle2)


# pro-pro

F_EV_PRO_PRO = mm.CustomNonbondedForce("eps*(sig/r)^6; eps=sqrt(eps1*eps2); sig = sig1 + sig2")
F_EV_PRO_PRO.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

F_EV_PRO_PRO.addPerParticleParameter('eps')
F_EV_PRO_PRO.addPerParticleParameter('sig')

for i in range (0, len(pro_pro_ex_vol)):
    eps = float(pro_pro_ex_vol[i][0])
    sig = float(pro_pro_ex_vol[i][1])

    F_EV_PRO_PRO.addParticle([eps*u.kilocalories_per_mole, sig*u.angstroms])
for j in range (0, len(exclusion_list)):
    particle1 = int(exclusion_list[j][0])
    particle2 = int(exclusion_list[j][1])

    F_EV_PRO_PRO.addExclusion(particle1, particle2)


# dna-pro



F_EV_DNA_PRO = mm.CustomNonbondedForce("delta(scale_fact)*eps*(sig/r)^6; scale_fact = scale_fact1 + scale_fact2; eps=sqrt(eps1*eps2); sig = sig1 + sig2")
F_EV_DNA_PRO.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

F_EV_DNA_PRO.addPerParticleParameter('eps')
F_EV_DNA_PRO.addPerParticleParameter('sig')
F_EV_DNA_PRO.addPerParticleParameter('scale_fact')

for i in range (0, len(dna_pro_ex_vol)):

    sig = float(dna_pro_ex_vol[i][0])
    eps = float(dna_pro_ex_vol[i][1])
    scale_fact = int(dna_pro_ex_vol[i][2])

    F_EV_DNA_PRO.addParticle([eps*u.kilocalories_per_mole, sig*u.angstroms, scale_fact])

for j in range (0, len(exclusion_list)):
    particle1 = int(exclusion_list[j][0])
    particle2 = int(exclusion_list[j][1])

    F_EV_DNA_PRO.addExclusion(particle1, particle2)



#~~~~~~~~~~~~~~~~~~~~~~~~Electrostatics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~DNA-DNA~~~~~~~~~#

F_Coulomb_DNA_DNA = mm.CustomNonbondedForce( "(coulomb_prefix/diele_dna)*charge*exp(-kappa_dna*r)/r; charge = charge1 * charge2;")
F_Coulomb_DNA_DNA.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

F_Coulomb_DNA_DNA.addPerParticleParameter('charge')
F_Coulomb_DNA_DNA.addGlobalParameter("coulomb_prefix", coulomb_prefix_val*u.kilocalories_per_mole*u.angstroms)
F_Coulomb_DNA_DNA.addGlobalParameter("kappa_dna", kappa_dna_val*u.angstroms**-1)
F_Coulomb_DNA_DNA.addGlobalParameter("diele_dna", diele_dna)

for i in range (0, len(dna_dna_DH)):
    charge = float(dna_dna_DH[i][0])
    F_Coulomb_DNA_DNA.addParticle([charge])

for j in range (0, len(exclusion_list)):
    particle1 = int(exclusion_list[j][0])
    particle2 = int(exclusion_list[j][1])

    F_Coulomb_DNA_DNA.addExclusion(particle1, particle2)


#~~~~~~~~PRO-PRO~~~~~~~~~#

F_Coulomb_PRO_PRO = mm.CustomNonbondedForce( "(coulomb_prefix/diele_pro)*charge*exp(-kappa_pro*r)/r; charge = charge1 * charge2;")
F_Coulomb_PRO_PRO.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

F_Coulomb_PRO_PRO.addPerParticleParameter('charge')
F_Coulomb_PRO_PRO.addGlobalParameter("coulomb_prefix", coulomb_prefix_val*u.kilocalories_per_mole*u.angstroms)
F_Coulomb_PRO_PRO.addGlobalParameter("kappa_pro", kappa_pro_val*u.angstroms**-1)
F_Coulomb_PRO_PRO.addGlobalParameter("diele_pro", diele_pro)

for i in range (0, len(pro_pro_DH)):
    charge = float(pro_pro_DH[i][0])
    F_Coulomb_PRO_PRO.addParticle([charge])

for j in range (0, len(exclusion_list)):
    particle1 = int(exclusion_list[j][0])
    particle2 = int(exclusion_list[j][1])

    F_Coulomb_PRO_PRO.addExclusion(particle1, particle2)



#~~~~~~~~DNA-PRO~~~~~~~~#

F_Coulomb_DNA_PRO = mm.CustomNonbondedForce( "delta(scale)*(coulomb_prefix/diele_dna_pro)*charge*exp(-kappa_dna_pro*r)/r; charge = charge1 * charge2; scale = scale1 + scale2")
F_Coulomb_DNA_PRO.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

F_Coulomb_DNA_PRO.addPerParticleParameter('charge')
F_Coulomb_DNA_PRO.addPerParticleParameter('scale')
F_Coulomb_DNA_PRO.addGlobalParameter("diele_dna_pro", diele_dna_pro)
F_Coulomb_DNA_PRO.addGlobalParameter("coulomb_prefix", coulomb_prefix_val*u.kilocalories_per_mole*u.angstroms)
F_Coulomb_DNA_PRO.addGlobalParameter("kappa_dna_pro", kappa_dna_pro_val*u.angstroms**-1)

for i in range (0, len(dna_pro_DH)):
    charge = float(dna_pro_DH[i][0])
    scale = int(dna_pro_DH[i][1])
    F_Coulomb_DNA_PRO.addParticle([charge, scale])
for j in range (0, len(exclusion_list)):
    particle1 = int(exclusion_list[j][0])
    particle2 = int(exclusion_list[j][1])

    F_Coulomb_DNA_PRO.addExclusion(particle1, particle2)

#~~~~~~~~~~~~~Add all the forces here~~~~~~~~~~~~~~#

system.addForce(F_Bond_DNA)
system.addForce(F_Bond_PRO)
system.addForce(F_Angle_DNA)
system.addForce(F_Angle_PRO)
system.addForce(Stackingforce)
system.addForce(Hbondingforce)
system.addForce(F_Native_PRO)
system.addForce(F_Native_DNA_PRO)
system.addForce(F_EV_DNA_DNA)
system.addForce(F_EV_PRO_PRO)
system.addForce(F_EV_DNA_PRO)
system.addForce(F_Coulomb_DNA_DNA)
system.addForce(F_Coulomb_PRO_PRO)
system.addForce(F_Coulomb_DNA_PRO)
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Simulation Block~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############~~~~~~~Group Forces~~~~~~~~~#################

for i in range(system.getNumForces()):
    force = system.getForce(i)
    force.setForceGroup(i)

###################~~~~~~~~Pick the integrator~~~~~~~~###################

integrator = mm.LangevinIntegrator(temp*u.kelvin, gamma_sim/u.picosecond, dt*u.femtosecond)


###################~~~~~~~~Simulation block~~~~~~~~###################

platform = mm.Platform.getPlatformByName(platform_name)
if platform_name == "CUDA":
    properties = {'Precision': 'double'}
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
else:
    simulation = app.Simulation(pdb.topology, system, integrator, platform)

simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temp*u.kelvin)

###################~~~~~~~~If restart is applicable then read the position and velocities~~~~~~~~###################


if restart==True:
    with open(prev_restart_filename, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())


######################Print Energies###################


print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

###################~~~~~~~~Report the simulation progress~~~~~~~~###################


datafile = open(datafile_name, 'w')
simulation.reporters.append(app.StateDataReporter(datafile, data_interval, step=True, potentialEnergy=True,  temperature=True))
simulation.reporters.append(EnergyCompReporter(energy_filename, data_interval))
simulation.reporters.append(app.DCDReporter(dcd_filename, data_interval_dcd))
simulation.reporters.append(app.CheckpointReporter(next_restart_filename, data_interval_dcd))
simulation.step(numsteps)


###################~~~~~~~~Report the time taken by the simulation~~~~~~~~###################

time_en = datetime.datetime.utcnow()
time_lapse = time_en - time_st
print('Simulation took\t')
print(time_lapse)
