#!/usr/bin/env python

# Import modules #

from PE5_header import *
from PE5_protein import Protein

########################
### Population class ###
########################

class Population:
	def __init__(self):
		self.popsize = None
		self.ParentGeneration = []
		self.CurrentGeneration = []
		self.NextGeneration = []
		self.CumFitness = []

	def Initialize(self):
		self.popsize = population_size
		self.CumFitness = np.zeros(self.popsize)
		for i in range(self.popsize):
			self.CurrentGeneration.append(Protein(i+1))
			self.CurrentGeneration[i].parent_idx = 0
			if genotype_level == 'nt':
				self.CurrentGeneration[i].nt_sequence = input_sequence
				self.CurrentGeneration[i].Translate()
			elif genotype_level == 'aa':
				self.CurrentGeneration[i].aa_sequence = input_sequence
			self.CurrentGeneration[i].MakePhenotype()
			if comp_scale == 0:   #Global competition.
				self.CurrentGeneration[i].AssignFitness(target_structure)

	def CalculateFitness(self):
		m = min([self.CurrentGeneration[i].fitness for i in range(self.popsize)])	#Since fitness is relative, we subtract the lowest fitness from all others, in order to not overflow on the exponentiation of e.
		for i in range(self.popsize):
			self.CumFitness[i] = (self.CumFitness[i-1] if i>0 else 0) + math.e ** (selection_coefficient * (self.CurrentGeneration[i].fitness-m))
		for i in range(self.popsize):
			self.CumFitness[i] = self.CumFitness[i] / self.CumFitness[-1]

	def GlobalCompetition(self):
		throw = rn.random()
		c = 0
		while c < len(self.CumFitness) and self.CumFitness[c] < throw:
			c += 1

		return c

	def LocalCompetition(self):
		candidates = rn.sample(range(self.popsize), k=comp_scale)

		if rn.random() < p_select_for_target:
			fits = [self.CurrentGeneration[c].SimilarityToPhenotype(target_structure) for c in candidates]
		else:
			fits = [self.CurrentGeneration[c].SimilarityToPhenotype(self.ParentGeneration[self.CurrentGeneration[c].parent_idx]) for c in candidates]

		max_f = [c for c, f in zip(candidates, fits) if f==max(fits)]
		return rn.choice(max_f)

	def Update(self):

		if comp_scale == 0:
			self.CalculateFitness()

		for i in range(self.popsize):
			if comp_scale == 0:
				c = self.GlobalCompetition()
			else:
				c = self.LocalCompetition()

			self.NextGeneration.append(Protein(i+1))
			self.NextGeneration[i].Replicate(self.CurrentGeneration[c])
			self.NextGeneration[i].Mutate()
			self.NextGeneration[i].MakePhenotype()
			if comp_scale == 0:
				self.NextGeneration[i].AssignFitness(target_structure)

		self.ParentGeneration = self.CurrentGeneration
		self.CurrentGeneration = self.NextGeneration
		self.NextGeneration = []

	def Output(self, time):
		with open(output_file, 'a') as fout:
			for i in range(self.popsize):
				fout.write(f"{time}\t{self.CurrentGeneration[i].PrintString()}\n")