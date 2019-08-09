#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
#import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import pickle
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import copy
from scipy.stats import sem, t
from scipy import mean
import re
import os
from time import sleep

#VARIABLES
Combinations= ["CancerImmune", "MetabolicImmune", "CancerMetabolic"]
Pairs_of_diseases = [["Cancer", "Immune"], ["Metabolic", "Immune"], ["Cancer", "Metabolic"]]
Pleiotropies=["Agon_early_early", "Agon_early_late", "Agon_early_early",
"Antagon_early_early", "Antagon_early_late", "Agon_early_early"]

#LOOPS
for age in range(10, 61):
	#First of all, get frequency and multiply to know probability of AP in that gene
	for pair in Pairs_of_diseases:
		INPUT_PATH1="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/" + pair[0] + "/results/Onset/"
		INPUT_PATH2="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/" + pair[1] + "/results/Onset/"
		get_frequency_of_AP_EarLate(INPUT_PATH1, age, pair[0])


def get_frequency_of_AP_EarLate(PATH, age, pair_value):
	Disease_path = PATH + pair_value + "_snps/Age_threeshold_" + str(age) + "/Antagon_early_late"
	with open(Disease_path, "r") as f:
   		size=sum(1 for _ in f)
		print(size)







	for filename in os.listdir(PATH)
