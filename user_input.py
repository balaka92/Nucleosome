# Simulation:

temp = 298.0                    #Temperature in Kelvin
salt_conc = 100                 #monovalent salt conc in mM
dt = 2.0                        #time step of simulation in fs
gamma_sim = 0.1                 #Friction in picoseconds**-1
numsteps = 20000                    #Number of steps of the simulation
numsteps_minimize = 0            #Number of steps for minimization
data_interval = 10       #Frequency at which energy components are saved
data_interval_dcd = 10   #Frequency at wich dcd files are saved


# restarting the simulation

restart = False          #True if you want to restart

#  input filenames:
pdb_filename = "INIT.pdb"  #PDB file consisting of DNA and Protein chains in a specific format, look at example.pdb
                           #DNA sequence has chain ID A,B and Histone proteins from C to J#
                           #Maintain the same scheme#

# output filenames:

datafile_name = "output_run_0.dat"                   #Change filename for every restart
prev_restart_filename = "chk0.chk"           #Change filename for every restart
next_restart_filename = "chk1.chk"           #Change filename for every restart
energy_filename = "energy_components_run_0.txt"     #Change filename for every restart
dcd_filename = "run_0.dcd"  #Change filename for every restart

#  platform:

platform_name = "CPU"              #options are "CPU" or "CUDA"
