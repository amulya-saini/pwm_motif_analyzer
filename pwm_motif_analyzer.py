#!/usr/bin/env python
# coding: utf-8

# In[1]:


# ignoring warnings
import warnings
warnings.filterwarnings('ignore')

# importing pandas package as pd
import pandas as pd

# importing Numpy package as np
import numpy as np


# In[2]:


# importing Library
import os

# printing the current directory
datapath = os.getcwd()


# In[3]:


# Joining the required file to the datapath
file1_path = datapath + "\\E_coli_K12_MG1655.400_50"

# print
# file1_path

# Loading the ppt file as pandas dataframe using read_csv function
df = pd.read_csv(file1_path, sep='\\', header=None)

# dropping the column with "NaN" values
df = df.drop(2, axis = 1)

# renaming the columns into gene_id and sequence
df = df.rename(columns = {0:"gene_id",1:"sequence"})

# adding a column "sequence_length" to calculate the length of sequences
df['sequence_length'] = df['sequence'].apply(len)

# checking if the sequence length is 452 and removing the space before and after sequence
if (df['sequence_length'] == 452).any():
    df['sequence'] = df['sequence'].str.strip()
    
# calculating the length again to make sure space is removed
df['sequence_length'] = df['sequence'].apply(len)

# removing space before or after gene_id
df['gene_id'] = df['gene_id'].str.strip()

# print
df


# In[4]:


# Joining the required file to the datapath
file2_path = datapath + "\\argR-counts-matrix.txt"

# print
# file2_path

# Loading the ppt file as pandas dataframe using read_csv function
df1 = pd.read_csv(file2_path, sep = "\t", header = None)

# dropping the column with values "|"
df1 = df1.drop(1, axis = 1)

# setting a,c,g,t as index
df1 = df1.set_index(0)

# converting the values from string to numeric to add +1 
df1 = df1.apply(pd.to_numeric, errors='coerce')

# adding +1 to all values to calculate log later
df1 = df1.add(1)

# Calculting position frequency matrix
    # by diving the each value with sum of the column
    # limiting the decimal to 3
pfm = df1.div(df1.sum()).round(3)

# print("pfm:", pfm)

# calculating position weight matrix by taking formula
    # "log2(frequency matrix value/0.25)"
pwm = np.log2(pfm/0.25)

# print("pwm:", pwm)

# calculating the sliding window size
window_size = pwm.shape[1]

# coverting the dataframe into a numpy array
pwm = pwm.to_numpy()

# Transposing the array to later correlate the index with base index in a sequence
pwm_t = np.transpose(pwm)

# print
pwm_t


# In[5]:


# creating a list of nucleotides present in the sequence
nucleotides = ["a","c","g","t"]


# In[6]:


# creating a function to calculate pwm score of a subsequence from a sequence
def similarity_score(subseq, pwm_t):
    
    # Creating a list of subsequence index looping through each nucleotide in the subsequence
    subseq_idx = [nucleotides.index(nucleotide) 
                  for nucleotide in subseq]
    # print(subseq_idx)
    
    # Summing the values from transposed pwm obtained using the length and index of subsequence
    return np.sum(pwm_t[np.arange(len(subseq)), subseq_idx])


# In[7]:


# creating an empty list to store the gene_ids and their scores
top_scores = []

# creating a for loop to iterate through every row in dataframe
for i, row in df.iterrows():
    
    # assigning variables to every column in the dataframe for every iteration
    gene_id = row['gene_id']
    sequence = row['sequence']
    sequence_length = row['sequence_length']
    
    # assigning "-infinity" as maximum score as the score can sum upto negative values
    max_score = -np.inf
    
    # creating a for loop to iterate through every subsequence from a given sequence
    for seq in range(sequence_length - window_size + 1):
        
        # defining subsequence by using ":" and predefined window size(18)
        subseq = sequence[seq : seq + window_size]
        
        # using the predefined similarity score function to calculate the score using pwm
        score = similarity_score(subseq, pwm_t)
        
        # if loop to check if the score is more than predefined max_score = -infinity
        if score > max_score:
            
            # replacing the max_score with the new max_score
            max_score = score
    
    # creating an if loop to check if the top scores list contains highest max scores
    # checking if the list contains less than 30 values
    if len(top_scores) < 30:
        
        # if the list is less than 30, adding the gene_id and max_score to the list 
        top_scores.append((gene_id, max_score))
            
    else:
        
        # creating a new variable that has the first value from top_scores list
        min_score_index = 0
        
        # the for loop then compares each score in the list to the score of min_score_index
        for a in range(1, 30):
            
            # the if loop checks the top_score and if it is lower than min_score index
            if top_scores[a][1] < top_scores[min_score_index][1]:
                
                # if the score is lower, it replaces the min_score_index with the lower value and it's gene id
                min_score_index = a
                
        # if loop to check if a gene id's max score lower than other gene's max score        
        if max_score > top_scores[min_score_index][1]:
            
            # replacing it with the highest score and it's gene_id
            top_scores[min_score_index] = (gene_id, max_score)
            
            
# sorting the list using lamda function to consider only the second element in the tuple (max score)
top_scores.sort(key = lambda x : x[1], reverse=True)

# getting top 30 gene ids using indexing after sorting
top30_gene_ids = [x[0] for x in top_scores]

# printing the top 30 ids
print(top30_gene_ids)


# In[8]:


# checking the length of the obtained gene_ids
len(top30_gene_ids)

