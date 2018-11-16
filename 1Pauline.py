# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:28:49 2018

@author: Yang Pauline Lilia
"""

import numpy as np

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
        print (self.__dict__)


'''
-------Input data----------
'''
data  = open("tousgenesidentiques.gff", 'r')
lines = data.readlines()
genes = []
for l in lines:
    genes.append(l.split())
data.close()
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
        a.orientation = 0
    gene_list.append(a)


# check genome
#for g in gene_list:
    #g.display()


#----------------------------


                


class Genome:
    '''Genome contains a set of genes, build in function is to simulate the 'evolution'
        
        Attributes:
            gene_list (list): a list contain all the gene.
            rate (float): the prob. ratio between insertion/deletion and inversion.
            generation (int): number of generations.
            genome_len(int): in the first generation length=30000
            prot_posi(int list): list of protein position,with defaut value
    '''
    def __init__(self):
        # each gene is a list [gene_id, start position, end positon, orientation]
        self.gene_list = []
        self.change_rate = 0.5
        self.generation = 0
        # defaut value is according to the TP 
        self.genome_len = 30000
        self.prot_posi = [1,3001,6001,9001,12001,15001,18001,21001,24001,27001]
    
    def insert(self, position=None,insert_len=60):
        '''insert a sequence in the genome
       
           Arguments:
               position (int/int list): insert position, defaut random choose
               insert_len (int): the length of insertion, defaut=60
       '''
       # if  there is no input position, random choose one 
        if position == None:
           r_position = self.select_modify_position()
        else:
           r_position = position
        #print "insert position: "+str(r_position)
        # find all genes after r_position, every gene + insert_len
        genes_to_modify = [g for g in self.gene_list if g.start > r_position]
        for g in genes_to_modify:
            g.start = g.start + insert_len
            g.end = g.end + insert_len
        
        #genome_len + insert_len
        self.genome_len += insert_len
    
    def delete(self, position=None, delete_len=60):
        '''delete a sequence in the genome
       
           Arguments:
               position (int): delete position, defaut random choose
               delete_len (int): the length of deletion, defaut=60
       '''
        if position == None:
           r_position = self.select_modify_position() #!!! not rubust need to change
           #sition: "+str(r_position)
           
        else:
           r_position = position
        
        # find all genes after r_position, every gene - delete_len
        genes_to_modify = [g for g in self.gene_list if g.start > r_position]   
        for g in genes_to_modify:
            g.start = g.start - delete_len
            g.end = g.end - delete_len
        
        self.genome_len -= delete_len

    def inversion(self, position1= None, position2=None): 
        Number_gene=0
        start_gene=[]
        end_gene=[]
        name_gene=[]
        orientation=[]
        if position1 == None:
            r_position1 = self.select_modify_position() #!!! not rubust need to change
            #print "delete position: "+str(r_position1)   
        else:
            r_position1 = position1
        if position2 == None:
            r_position2 = self.select_modify_position() #!!! not rubust need to change
            #print "delete position: "+str(r_position2)   
        else:
            r_position2 = position2   
        # inversion between r_posi1 and r_posi2
        for g in self.gene_list:
            if (r_position1 <= g.start and r_position2>= g.end) : 
                Number_gene = Number_gene +1
                start_gene.append(g.start)
                end_gene.append(g.end)
                name_gene.append(g.id)
                orientation.append(g.orientation)
        start_gene_i=min(start_gene)
        end_gene_i=max(end_gene)
        dist_gene=[]
        n=len(name_gene)#nombre de gene 
        l=1000 # taille d'un gene 
        for i in range(1,n):
        		dist_gene.append(start_gene[i]-end_gene[i-1])
        	#inverser les index des nouveaux genes 
        name_gene = name_gene[::-1]
        dist_gene =dist_gene[::-1]
        	#ecrire la nouvelle premiere position (start et end) 
        start_gene[0]= r_position1 + r_position2- end_gene_i
        end_gene[0]= start_gene[0] + l
        orientation[0]=  orientation[0]*(-1)
        for i in range(1,n): 
        		start_gene[i]= end_gene[i-1] + dist_gene[i-1]
        		end_gene[i]= start_gene[i]+l 
        		orientation[i]=  orientation[i]*(-1)
        	
            #g_to_modify=[gn for gn in self.gene_list if gn.id in name_gene]
        gene_to_modify = list(zip(name_gene, start_gene, end_gene,orientation))
        gene_to_modify.sort(key= lambda g:g[0])
        #print (gene_to_modify)
        
        for i,g in enumerate([gn for gn in self.gene_list if gn.id in name_gene]):
            		g.start = gene_to_modify[i][1]
        		#print g.start
            		g.end = gene_to_modify[i][2]
           		g.orientation = gene_to_modify[i][3]
        		#print name_gene[i], g.id
        
        	#print self.gene_list		
        return name_gene, start_gene, end_gene , orientation
	
	
	
    def select_modify_position(self):
        ''' randomly select a position can be modify
            
            Return a position in the Genome is modifible
        '''
        # identify the positions we can not touch
        untouchble = []
        # nothing can be insert inside a gene
        for g in self.gene_list:
            untouchble+=range(g.start,g.start+g.length*g.orientation,g.orientation)
            
        genome_all_posi = range(1,1+self.genome_len)
        modifible = list(set(genome_all_posi)-set(untouchble))
        r_position = np.random.choice(modifible)                
            
        return r_position
        #print untouchble

gn1 = Genome()
gn1.gene_list=gene_list
#print gn1.gene_list
#gn1.insert(insert_len=6000)
#gn1.delete(delete_len=3000)
[name, start, end, prientation]=gn1.inversion(position1=5, position2=24000)
#print name, start
#gn1.update_inversion()
for g in gn1.gene_list:
	g.display()

  