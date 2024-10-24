#!/usr/bin/env python

# Define command line-adjustable parameters.

import Header

project_name = 'test'
initial_seed = 203
population_size = 100
simulation_time = 100
input_sequence = ''

# Variables (part 2) #

target_structure = ""
selection_coefficient = 0.1
comp_scale = 0
p_select_for_target = 0.1
fitness_criterion = 'target_structure'

# Variables (part 3) #

mutation_rate = 0.001
p_insertion = 0.
p_deletion = 0.
genotype_level = 'nt'
phenotype_level = '2D'
exp_threshold = 0.5