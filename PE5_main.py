#!/usr/bin/env python

# (24/10/2024)	Simulate protein evolution with various fitness options, but now using object-oriented programming.

#############
### Setup ###
#############

from PE5_header import *
from PE5_protein import Protein
from PE5_population import Population

def PrintHelp():
	print("-> Incorrent programme call\n::: Protein Evolution v.5 :::\nUsage: ./PE5_main.py -p [project name] -s [initial seed] -m [mutation rate] -w [selection coefficient] -N [population size] -t [simulation time] -F [fitness criterion] -C [competition scale] -T [target structure file] -i [input sequence]\nSee PE5_header.py for default parameters and other important settings...\n")
	sys.exit(1)

i_arg = 1
while i_arg != len(sys.argv):
	print(i_arg, sys.argv[i_arg])
	if sys.argv[i_arg] == '-p':
		print(project_name)
		change_project_name(sys.argv[i_arg+1])
		print(project_name)
	elif sys.argv[i_arg] == '-s':
		print(initial_seed)
		change_initial_seed(int(sys.argv[i_arg+1]))
		print(initial_seed)
	elif sys.argv[i_arg] == '-t':
		simulation_time = int(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-m':
		change_mutation_rate(float(sys.argv[i_arg+1]))
	elif sys.argv[i_arg] == '-w':
		change_selection_coefficient(float(sys.argv[i_arg+1]))
	elif sys.argv[i_arg] == '-N':
		change_population_size(int(sys.argv[i_arg+1]))
	elif sys.argv[i_arg] == '-F':
		change_fitness_criterion(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-C':
		change_comp_scale(int(sys.argv[i_arg+1]))
	elif sys.argv[i_arg] == '-T':
		change_target_structure(sys.argv[i_arg+1])
	elif sys.argv[i_arg] == '-i':
		change_input_sequence(sys.argv[i_arg+1])
	else:
		PrintHelp()
	i_arg += 2

if project_name == "" or initial_seed == 0:
	PrintHelp()
else:
	rn.seed(initial_seed)
	proj_dir = '/linuxhome/tmp/sam/protevol/' + project_name
	os.system(f'mkdir -p {proj_dir}')
	os.system(f'mkdir -p {proj_dir}/tmp')
	change_output_file(f'{proj_dir}/{project_name}.out')

if fitness_criterion == 'target_structure' and not target_structure:
	PrintHelp()

rn.seed(seed)

##################
### Simulation ###
##################

P = Population()
P.Initialize()
P.Output(0)

for time in range(1, simulation_time):
	P.Update()
	P.Output(time)