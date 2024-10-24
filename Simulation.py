#!/usr/bin/env python

# (24/10/2024)	Simulate protein evolution with various fitness options, but now using object-oriented programming.

#############
### Setup ###
#############

from Header import *
import Config
from Protein import Protein
from Population import Population

def PrintHelp():
	print("-> Incorrent programme call\n::: Protein Evolution v.5 :::\nUsage: ./PE5_main.py -p [project name] -s [initial seed] -m [mutation rate] -w [selection coefficient] -N [population size] -t [simulation time] -F [fitness criterion] -C [competition scale] -T [target structure file] -i [input sequence] [-pst proportion selecting for target]\nSee PE5_header.py for default parameters and other important settings...\n")
	sys.exit(1)

i_arg = 1
while i_arg != len(sys.argv):
	if sys.argv[i_arg] == '-p':
		Config.project_name = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-s':
		Config.initial_seed = int(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-t':
		simulation_time = int(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-m':
		Config.mutation_rate = float(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-w':
		Config.selection_coefficient = float(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-N':
		Config.population_size = int(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-F':
		Config.fitness_criterion = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-C':
		Config.comp_scale = int(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-T':
		Config.target_structure = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-i':
		Config.input_sequence = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-pst':
		Config.p_select_for_target = float(sys.argv[i_arg+1])
	else:
		PrintHelp()
	i_arg += 2

if Config.project_name == "" or Config.initial_seed == 0:
	PrintHelp()
else:
	rn.seed(Config.initial_seed)
	os.system(f'mkdir -p {linuxhome_dir}/{Config.project_name}')
	os.system(f'mkdir -p {linuxhome_dir}/{Config.project_name}/tmp')
	os.system(f'touch {linuxhome_dir}/{Config.project_name}/{Config.project_name}')

if Config.fitness_criterion == 'target_structure' and not Config.target_structure:
	PrintHelp()

##################
### Simulation ###
##################

P = Population()
P.Initialize()
P.Output(0)

for time in range(1, simulation_time):
	print(time)
	P.Update()
	P.Output(time)