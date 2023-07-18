# pwm_motif_analyzer
Python-based tool for analyzing sequence motifs and identifying high-scoring subsequences based on position weight matrices (PWMs). This repository provides a collection of scripts and functions to process and analyze gene expression data and DNA sequences, with a focus on calculating PWM scores and identifying top-scoring subsequences.


Programming Language and Version: Python 3.9.13

Required packages/Libraries: 
1, Pandas
- To install Pandas package, type the following command in command prompt
	py -m pip install "pandas"
- To import the package to start using it, use the following command
	import pandas as pd

2, Numpy
- To install Numpy package, type the following command in command prompt
	py -m pip install "numpy"
- To import the package to start using it, use the following command
	import numpy as np

Required input files:
1, argR-counts-matrix.txt 
2, E_coli_K12_MG1655.400_50
 
Description: 
- The read_csv functions reads both the files respective seperators as a pandas data frame.
- The pandas and numpy package was used for further data processing.
- The user defined similarity_score function calculates the pwm score for each subsequence in sequence.
- The code then stores highest score of each gene id, then prints the top 30 gene ids after sorting the gene ids based on the their pwm scores. 

Execution: 
1,Open Command Prompt.
2, Change the working directory to the location of the file.
	cd path\to\the\file\location
3, Run the script using the following command.
	pwm_motif_analyzer.py
4, Once execution is complete check output.
	
Output Files:
None

Author: Amulya Saini
Date: 04/10/2023
