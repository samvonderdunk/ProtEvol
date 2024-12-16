#!/usr/bin/env python

# Import modules #

from Header import *
import Config
from Protein import Protein

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
		self.popsize = Config.population_size
		self.CumFitness = np.zeros(self.popsize)
		for i in range(self.popsize):

			if i == 0:	#We actually calculate all properties for the first individual.
				self.CurrentGeneration.append(Protein(i+1))
				self.CurrentGeneration[i].parent_idx = 0

				if Config.genotype_level == 'nt':
					self.CurrentGeneration[i].nt_sequence = Config.input_sequence
					self.CurrentGeneration[i].Translate()
				elif Config.genotype_level == 'aa':
					self.CurrentGeneration[i].aa_sequence = Config.input_sequence

				self.CurrentGeneration[i].MakePhenotype()

				if Config.comp_scale == 0:   #Global competition.
					self.CurrentGeneration[i].AssignFitness()

			else:	#For the remaining individuals, we just copy from the first individual.
				self.CurrentGeneration.append(Protein(i+1))
				self.CurrentGeneration[i].Replicate(self.CurrentGeneration[0])
				self.CurrentGeneration[i].parent_idx = 0
				self.CurrentGeneration[i].CopyPhenotype(self.CurrentGeneration[0])

	def CalculateFitness(self):
		m = min([self.CurrentGeneration[i].fitness for i in range(self.popsize)])	#Since fitness is relative, we subtract the lowest fitness from all others, in order to not overflow on the exponentiation of e.
		for i in range(self.popsize):
			self.CumFitness[i] = (self.CumFitness[i-1] if i>0 else 0) + math.e ** (Config.selection_coefficient * (self.CurrentGeneration[i].fitness-m))
		for i in range(self.popsize):
			self.CumFitness[i] = self.CumFitness[i] / self.CumFitness[-1]

	def GlobalCompetition(self):
		throw = rn.random()
		c = 0
		while c < len(self.CumFitness) and self.CumFitness[c] < throw:
			c += 1

		return c

	def LocalCompetition(self):
		candidates = rn.sample(range(self.popsize), k=Config.comp_scale)

		if rn.random() < Config.p_select_for_target:
			fits = [self.CurrentGeneration[c].SimilarityToPhenotype(Config.target_structure) - len(self.CurrentGeneration[c].aa_sequence)/3 for c in candidates]
		else:
			if len(self.ParentGeneration) == 0:	#First timestep, so there is no parent generation yet.
				fits = [1 for c in candidates]
			else:
				fits = [self.CurrentGeneration[c].SimilarityToPhenotype(self.ParentGeneration[self.CurrentGeneration[c].parent_idx-1].ss_structure) / len(self.CurrentGeneration[c].ss_structure) for c in candidates] #Parent similarity normalized by length of current individual to not incentivize sequence growth.

		max_f = [c for c, f in zip(candidates, fits) if f==max(fits)]
		return rn.choice(max_f)

	def Update(self):

		if Config.comp_scale == 0:
			self.CalculateFitness()

		for i in range(self.popsize):
			if Config.comp_scale == 0:
				c = self.GlobalCompetition()
			else:
				c = self.LocalCompetition()

			self.NextGeneration.append(Protein(i+1))
			self.NextGeneration[i].Replicate(self.CurrentGeneration[c])
			if self.NextGeneration[i].Mutate():	#No mutation happened, so just copy the phenotype
				self.NextGeneration[i].CopyPhenotype(self.CurrentGeneration[c])
			else:
				self.NextGeneration[i].MakePhenotype()
				if Config.comp_scale == 0:
					self.NextGeneration[i].AssignFitness()

		self.ParentGeneration = self.CurrentGeneration.copy()
		self.CurrentGeneration = self.NextGeneration.copy()
		self.NextGeneration = []

	def Output(self, time):
		with open(f'{linuxhome_dir}/{Config.project_name}/{Config.project_name}.out', 'a') as fout:
			for i in range(self.popsize):
				fout.write(f"{time}\t{self.CurrentGeneration[i].PrintString()}\n")