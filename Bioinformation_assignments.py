#!/usr/bin/env python
# coding: utf-8

# In[3]:


def k_mers(neucleotide_number):
    neucleotides = "AGTC"
    k_mers = ["A", "G", "T", "C"]
    k_mers_loop = []
    while len(k_mers[0]) <= neucleotide_number - 1:
        for i1 in neucleotides:
            for i2 in k_mers:
                k_mers_loop.append(i1+i2)
        k_mers = k_mers_loop
        k_mers_loop = []
    return k_mers


# In[4]:


def hamming_distance(a, b):
    d = 0
    for i,j in zip(a, b):
        if i != j:
            d += 1
    return d
            


# In[5]:


def DNA_strings_k_mers(DNA_strings, k_mers_length):
    DNA_strings_k_mers = []
    for j1 in range(0,len(DNA_strings)):
        DNA_strings_k_mers.append([])
        for j2 in range(0, len(DNA_strings[j1])-(k_mers_length)):
            j2_motif = DNA_strings[j1][j2 : j2 + k_mers_length]
            DNA_strings_k_mers[j1].append(j2_motif)
    return DNA_strings_k_mers


# In[6]:


def hamming_distance_sheet(DNA_strings_k_mers, k_mers):
    hamming_distance_sheet = []
    for k1 in range(0, len(k_mers)):
        hamming_distance_sheet.append([])
        for k2 in DNA_strings_k_mers:
            d = []
            for k3 in k2:
                d.append(hamming_distance(k_mers[k1], k3))
            hamming_distance_sheet[k1].append(min(d))
    return hamming_distance_sheet
                 


# In[7]:


def motif_finding(DNA_strings, k_mer_length):
    total_d_sum = []
    for l1 in range(0, len(hamming_distance_sheet(DNA_strings_k_mers(DNA_strings, k_mer_length), k_mers(k_mer_length)))):
         total_d_sum.append(sum(hamming_distance_sheet(DNA_strings_k_mers(DNA_strings, k_mer_length), k_mers(k_mer_length))[l1]))
    return k_mers(k_mer_length)[total_d_sum.index(min(total_d_sum))]

        


# In[ ]:


DNA_strings = ["AGGCTTCTGATCGCCGTAGGCTACGAGATCTTTAGCCTCCAC",
"AAATTTCGCGGGTTTAGCTGGGCGGTTCGACAGGATTCCTGT",
"TTTAGCCACCTGACGCGGGGACCCGGCTGGCCCTGAAACTTA",
"TTTAGTCCGTGCTAGCAAGTACTTCGAGTACCGCATTCAAAT",
"TGACTCTACGTAATAGGATTTAGTGAGGCCCGCCTGGCACCG",
"GATGTTTTTAGTACTGTACGAACCAGGTGTATTGCTTTGCGG",
"GAGCAGCCGAGGTGGGCTGATTATACAATCTTTAGAGGAGCA",
"TTTAGCAGTATTTGTCCCCCTGGTGATCCTACATTCAGCCGA",
"TAAGCTTGCGAGACCCCAGATCAGTCTTGTTTTAGCACCACC",
"GGACGAGCGTGGTTTAGGCAGTATGACCGCAGACACACGTTA"]

k_mer_length = 5
motif_finding(DNA_strings, k_mer_length)

