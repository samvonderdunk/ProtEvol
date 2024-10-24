#!/usr/bin/env python

import random as rn
import numpy as np
import torch
import math, sys, os
import KC

model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")

import esm
os.environ['MKL_THREADING_LAYER'] = 'GNU'	#Solves an issue with numpy.

model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
# model.set_chunk_size(128)

# Simulation variables #

linuxhome_dir = '/linuxhome/tmp/sam/protevol'
project_name = 'test'
initial_seed = 203
population_size = 100
simulation_time = 100
output_file = ''
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

# Fixed variables #

Nucleotides = ['A','C','G','U']
MutProbsNT = [[0, 1/6, 4/6, 1/6], [1/6, 0, 1/6, 4/6], [4/6, 1/6, 0, 1/6], [1/6, 4/6, 1/6, 0]]
AminoAcids = ['A','G','I','L','P','V','F','W','Y','D','E','R','H','K','S','T','C','M','N','Q']
MutProbsAA = np.ones((20,20))

CodonTable = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

#Table from wikipedia but source is Tien et al. 2013 (empirical data).
MaxExposureTable = {
	'A': 121,
	'R': 265,
	'N': 187,
	'D': 187,
	'C': 140,
	'E': 214,
	'Q': 214,
	'G': 97,
	'H': 216,
	'I': 195,
	'L': 191,
	'K': 230,
	'M': 203,
	'F': 228,
	'P': 154,
	'S': 143,
	'T': 163,
	'W': 264,
	'Y': 255,
	'V': 165
}


## Change functions (allows you to change the parameter globally from PE5_main.py ##

def change_project_name(new_value):
    global project_name
    project_name = new_value

def change_initial_seed(new_value):
    global initial_seed
    initial_seed = new_value

def change_population_size(new_value):
    global population_size
    population_size = new_value

def change_output_file(new_value):
    global output_file
    output_file = new_value

def change_input_sequence(new_value):
    global input_sequence
    input_sequence = new_value

def change_target_structure(new_value):
    global target_structure
    target_structure = new_value

def change_selection_coefficient(new_value):
    global selection_coefficient
    selection_coefficient = new_value

def change_comp_scale(new_value):
    global comp_scale
    comp_scale = new_value

def change_p_select_for_target(new_value):
    global p_select_for_target
    p_select_for_target = new_value

def change_fitness_criterion(new_value):
    global fitness_criterion
    fitness_criterion = new_value

def change_mutation_rate(new_value):
    global mutation_rate
    mutation_rate = new_value

def change_p_insertion(new_value):
    global p_insertion
    p_insertion = new_value

def change_p_deletion(new_value):
    global p_deletion
    p_deletion = new_value

def change_genotype_level(new_value):
    global genotype_level
    genotype_level = new_value

def change_phenotype_level(new_value):
    global phenotype_level
    phenotype_level = new_value

def change_exp_threshold(new_value):
    global exp_threshold
    exp_threshold = new_value