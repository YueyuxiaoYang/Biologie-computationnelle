# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:28:49 2018

@author: Yang Pauline Lilia
"""

import numpy as np

'''
-------Input data----------
'''
data  = open("tousgenesidentiques.gff", 'r')
lines = data.readlines()
genes = []
for l in lines:
    genes.append(l.split())
    
genes = genes[5:-1]
# save genes in a list [class_gene, ....]
gene_list = [] # data structure for saving initial genome
for g in genes:
    a = gene()
    a.id = int(g[-1][12:])
    a.start = int(g[3])
    a.end = int(g[4])
    if  g[-3] == '+':
        a.orientation = 1
    else:
        a.orientation = -1
    gene_list.append(a)

# check genome
for g in gene_list:
    g.display()


#----------------------------

class gene:
    '''gene is a data structure for a gene
    
        Attributes:
            id (int): gene's id or name, to distinguish genes.
            start (int): start position.
            end (int): end position.
            orientation (bool):  '+':1, '-':-1.
            length (int): in this project length is a constant=1000.
    
    '''
    def __init__(self):
         # define defaut value
        self.id = 0
        self.start = -1 # start position should be between 1 and 30000
        self.end = 0 # end position should be [1:30000]
        self.orientation = 0 # '+':1, '-':0
        self.length = 1000 # length of the gene is 1000 by defaut,and it is a constant 
    def display(self):
        '''
            display all Attributes of a gene
        '''
        print self.__dict__
                


class Genome:
    '''Genome contains a set of genes, build in function is to simulate the 'evolution'
        
        Attributes:
            gene_list (list): a list contain all the gene.
            rate (float): the prob. ratio between insertion/deletion and inversion.
            generation (int): number of generations.
            genome_len(int): in the first generation length=30000
    '''
    def __init__(self):
        # each gene is a list [gene_id, start position, end positon, orientation]
        self.gene_list = []
        self.change_rate = 0.5
        self.generation = 0
        self.genome_len = 30000
    
    def insert(self, position=None,insert_len=1):
       '''insert a sequence in the genome
       
           Arguments:
               position (int/int list): insert position
               seq (?): sequence to insert
               insert_len (int): the length of insertion
       '''
       # if  there is no input position, random choose one 
       if position == None:
           pass
           #position = 
    def select_modify_position(self):
        ''' randomly select a position can be modify
            
            Return a position in the Genome is modifible
        '''
        # identify the positions we can not touch
        untouchble = []
        for g in self.gene_list:
            untouchble+=range(g.start,g.start+g.length*g.orientation,g.orientation)
            
        genome_all_posi = range(1,1+self.genome_len)
        modifible = list(set(genome_all_posi)-set(untouchble))
            
        return modifible
        #print untouchble

gn1 = Genome()
gn1.gene_list = gene_list[1]

           
    
    
    
    
    