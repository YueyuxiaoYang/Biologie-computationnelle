# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:28:49 2018

@author: Yang Pauline Lilia
"""

import numpy as np

'''
-------Input data----------
'''
data  = open("/home/kevin1024/Desktop/sam mayer/tousgenesidentiques.gff", 'r')
lines = data.readlines()
genes = []
for l in lines:
    genes.append(l.split())

#----------------------------



'''
-------define genome class
'''
class Genome:
    
    def __init__(self):
        # each gene is a list [gene_id, start position, end positon, orientation]
        self.gene_list = []
    