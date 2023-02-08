import numpy as np
import pandas as pd
import scipy
from plotnine import *


'''
Homework 1 problem 8 -- global alignment
use the simple scoring method of +1 for match and -1 for mismatch and indel
print the global alignment with one line per string with '-' characters for indels
'''
def global_alignment(sequence1, sequence2):
    match = 1
    mismatch_indel = -1
    #initialize "array" for scoring, could have used numpy array for array ops but alas
    arr = [ [0 for x in range(len(sequence1)+1)] for i in range(len(sequence2)+1)]
    
    #set horizontal and vertical col/row to -1 times index
    for index, num in enumerate(arr[0]): 
        arr[0][index] = mismatch_indel * index
    for index, num in enumerate(arr): 
        arr[index][0] = mismatch_indel * index

    #score array by adding +1 on the diagonal if chars are a match else -1 from the highest number around it
    for vert_index in range(1, len(arr)): 
        for horizontal_index in range(1, len(arr[vert_index])):
            if sequence1[horizontal_index-1] == sequence2[vert_index-1]:
                arr[vert_index][horizontal_index] = arr[vert_index-1][horizontal_index-1]+1
            else: 
                arr[vert_index][horizontal_index] = max(arr[vert_index-1][horizontal_index-1], 
                                                        arr[vert_index-1][horizontal_index], 
                                                        arr[vert_index][horizontal_index-1])-1
    #Backtrace 
    stopping = False #bool keeping track of when to stop
    pos = [len(sequence1), len(sequence2)] #keeping in mind seq1 is horizontal len and seq2 is vertical
    aligned_seq1 = ""
    aligned_seq2 = ""
    while not stopping: 
        box_max = max(arr[pos[1]-1][pos[0]-1], arr[pos[1]-1][pos[0]], arr[pos[1]][pos[0]-1])
        if box_max == arr[pos[1]-1][pos[0]-1]:
            aligned_seq1 = sequence1[pos[0]-1] + aligned_seq1
            aligned_seq2 = sequence2[pos[1]-1] + aligned_seq2
            #move diagonally
            pos[0] = pos[0]-1
            pos[1] = pos[1]-1
        #move horizontal
        elif box_max == arr[pos[1]][pos[0]-1]: 
            aligned_seq2 = "-" + aligned_seq2
            aligned_seq1 = sequence1[pos[0]-1] + aligned_seq1
            pos[0] = pos[0]-1
        #move vertical
        else:
            aligned_seq2 = sequence2[pos[1]-1] + aligned_seq2
            aligned_seq1 = "-" + aligned_seq1
            pos[1] = pos[1]-1
        #check whether or not to stop            
        if pos == [0,0]:
            stopping = True
    print("\n" + aligned_seq1)
    print(aligned_seq2)




'''
support code for creating random sequence, no need to edit
'''
def random_sequence(n):
    return("".join(np.random.choice(["A","C","G","T"], n)))

'''
support code for mutating a sequence, no need to edit
'''
def mutate(s, snp_rate, indel_rate):
    x = [c for c in s]
    i = 0
    while i < len(x):
        if np.random.random() < snp_rate:
            x[i] = random_sequence(1)
        if np.random.random() < indel_rate:
            length = np.random.geometric(0.5)
            if np.random.random() < 0.5: # insertion
                x[i] = x[i] + random_sequence(length)
            else:
                for j in range(i,i+length):
                    if j < len(x):
                        x[j] = ""
                    i += 1
        i += 1
    return("".join(x))

# creating related sequences
s1 = random_sequence(100)
s2 = mutate(s1, 0.1, 0.1)

# running your alignment code
global_alignment(s1, s2)
