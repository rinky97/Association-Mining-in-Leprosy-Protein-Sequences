# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 22:46:08 2018

@author: Rinky
"""

#importing libraries
import numpy as np
import pandas as pd
from Bio import SeqIO
import collections
from collections import Counter

#parsing the file
with open('sequence.fasta') as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    sequences = []
    A =[]
    C =[]
    D =[]
    E =[]
    F =[]
    G =[]
    H =[]
    I =[]
    K =[]
    L =[]
    M =[]
    N =[]
    P =[]
    Q =[]
    R =[]
    S =[]
    T =[]
    U =[]
    V =[]
    W =[]
    X =[]
    Y =[]
    a =[]
    for seq_record in SeqIO.parse("sequence.fasta","fasta"):
        identifiers.append(seq_record.id)
        lengths.append(len(seq_record.seq))
        sequences.append(seq_record.seq)
        A.append(seq_record.seq.count('A'))
        C.append(seq_record.seq.count('C'))
        D.append(seq_record.seq.count('D'))
        E.append(seq_record.seq.count('E'))
        F.append(seq_record.seq.count('F'))
        G.append(seq_record.seq.count('G'))
        H.append(seq_record.seq.count('H'))
        I.append(seq_record.seq.count('I'))
        K.append(seq_record.seq.count('K'))
        L.append(seq_record.seq.count('L'))
        M.append(seq_record.seq.count('M'))
        N.append(seq_record.seq.count('N'))
        P.append(seq_record.seq.count('P'))
        Q.append(seq_record.seq.count('Q'))
        R.append(seq_record.seq.count('R'))
        S.append(seq_record.seq.count('S'))
        T.append(seq_record.seq.count('T'))
        U.append(seq_record.seq.count('U'))
        V.append(seq_record.seq.count('V'))
        W.append(seq_record.seq.count('W'))
        X.append(seq_record.seq.count('X'))
        Y.append(seq_record.seq.count('Y'))
 
#creating a database       
bio_df = pd.DataFrame(np.column_stack([identifiers, lengths, sequences,A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y]), 
                               columns=['idTitle', 'lenTitle', 'seqTitle', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

#pruning those sequences that contain U and X by index values
bio_df=bio_df.drop(bio_df.index[[1012,1013,1671,1672,1673,1674,1675,1676,1677,1678,1679,1680,1681,1682,1683,1684,1685,1686,1687,1688,1690,1691,2631,3004,]])

#chars stores each sequence string as characters
x = bio_df['seqTitle']
chars = []
for line in x:
    chars.append(list(line))
       
#dictionary to store the occurance of each amino acid
amino_count ={}
for i in chars:
    for j in i:
        if j in amino_count:
            amino_count[j] += 1
        else:
            amino_count[j] = 1          
s = sum(bio_df['lenTitle'])
               
#dictionary that stores the normal average of each amino acid in a sequence
sup_amino={}
for key, values in amino_count.items():
    sup_amino[key] = values/s
    
freq_1_items = {}
for key, values in sup_amino.items():
    if(values >= 0.05):
        freq_1_items[key] = values
#a = bio_df[0]
#arr=[bio_df[3:]]


from itertools import combinations      #for finding out all possible aminoacid combinations
#arr = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
arr2=[]
for key,value in freq_1_items.items():
    arr2.append(key)
#finding the count of 2-amino acid pairs
comb_2 = list(combinations(arr2, 2))
pair2_count = []
for each in comb_2:
    pair2_count.append(list(map(min, zip(bio_df[each[0]], bio_df[each[1]]))))
pair2_sum = []
for each in pair2_count:
    pair2_sum.append(sum(each))
    

dict_2_amino = dict(zip(comb_2, pair2_sum))

sup_2_amino = {}
for key, values in dict_2_amino.items():
    sup_2_amino[key] = values/s

freq_2_items = {}
for key, values in sup_2_amino.items():
    if(values >= 0.05):
        freq_2_items[key] = values
    
#finding the count of 3-amino acid pairs
arr3 = []
for key,value in freq_2_items.items():
    for i in key:
        arr3.append(i)
arr3 = list(set(arr3))

comb_3 = list(combinations(arr3,3))
pair3_count = []
for each in comb_3:
    pair3_count.append(list(map(min, zip(bio_df[each[0]], bio_df[each[1]], bio_df[each[2]]))))
pair3_sum = []
for each in pair3_count:
    pair3_sum.append(sum(each))
    
dict_3_amino = dict(zip(comb_3, pair3_sum))

sup_3_amino = {}
for key, values in dict_3_amino.items():
    sup_3_amino[key] = values/s
    
freq_3_items = {}
for key, values in sup_3_amino.items():
    if(values >= 0.05):
        freq_3_items[key] = values
    
#finding the count of 4-amino acid pairs
arr4 = []
for key,value in freq_3_items.items():
    for i in key:
        arr4.append(i)
arr4 = list(set(arr4))        
        
comb_4 = list(combinations(arr4,4))
pair4_count = []
for each in comb_4:
    pair4_count.append(list(map(min, zip(bio_df[each[0]], bio_df[each[1]], bio_df[each[2]],bio_df[each[3]]))))
pair4_sum = []
for each in pair4_count:
    pair4_sum.append(sum(each))
    
dict_4_amino = dict(zip(comb_4, pair4_sum))

sup_4_amino = {}
for key, values in dict_4_amino.items():
    sup_4_amino[key] = values/s
    
freq_4_items = {}
for key, values in sup_4_amino.items():
    if(values >= 0.05):
        freq_4_items[key] = values
#finding the count fo 5-amino acid pairs    
        
arr5 = []
for key,value in freq_4_items.items():
    for i in key:
        arr5.append(i)
arr5 = list(set(arr5))

comb_5 = list(combinations(arr5,5))
pair5_count = []
for each in comb_5:
    pair5_count.append(list(map(min, zip(bio_df[each[0]], bio_df[each[1]], bio_df[each[2]], bio_df[each[3]], bio_df[each[4]]))))
pair5_sum = []
for each in pair5_count:
    pair5_sum.append(sum(each))
    
dict_5_amino = dict(zip(comb_5, pair5_sum))

sup_5_amino = {}
for key, values in dict_5_amino.items():
    sup_5_amino[key] = values/s
    
freq_5_items = {}
for key, values in sup_5_amino.items():
    if(values >= 0.05):
        freq_5_items[key] = values

