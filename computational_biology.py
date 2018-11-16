# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:28:49 2018

@author: Yang Pauline Lilia
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

''' if you wan to import this file 

sys.path
# add the path of this file in sys so as to import it
sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle') 
import computational_biology as cb
'''


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
        print (self.__dict__)
                


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
           r_position = self.select_modify_position_insert()
        else:
           r_position = position
        print ("insert position: "+str(r_position))
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
           r_position = self.select_modify_position_delete(delete_len) #!!! not rubust need to change
           print ("delete position: "+str(r_position))
           
        else:
           r_position = position
        
        # find all genes after r_position, every gene - delete_len
        genes_to_modify = [g for g in self.gene_list if g.start > r_position]   
        for g in genes_to_modify:
            g.start = g.start - delete_len
            g.end = g.end - delete_len
        
        self.genome_len -= delete_len
       
          
    def select_modify_position_insert(self):
        ''' randomly select a position can be modify
            
            Return a position in the Genome is modifible
        '''
        # identify the positions we can not touch
        untouchble = []
        # nothing can be insert inside a gene
        for g in self.gene_list:
            untouchble+=range(g.start,g.end+1)
            
        genome_all_posi = range(1,1+self.genome_len)
        modifible = list(set(genome_all_posi)-set(untouchble))
        r_position = np.random.choice(modifible)                
                            
        return r_position
        #print untouchble
    def select_modify_position_delete(self,delete_len):
        ''' randomly select a position can be modify
            
            Return a position in the Genome is modifible, 
            constrain:  
                delete position - delete_len not in the gene
        '''
        # identify the positions we can not touch
        untouchble = []
        # delete_posi - delete_len not in the gene
        for g in self.gene_list:
            untouchble+=range(g.start,g.end+delete_len+1)
        genome_all_posi = list(range(1+delete_len,1+self.genome_len))
        modifible = list(set(genome_all_posi)-set(untouchble))
        r_position = np.random.choice(modifible)   
        return r_position
    def select_modify_position_inversion(self):
        ''' randomly select a position can be modify
            
            Return 2 positions in the Genome is modifible
        '''
        # identify the positions we can not touch
        untouchble = []
        # delete_posi - delete_len not in the gene
        for g in self.gene_list:
            untouchble+=range(g.start,g.end+1)
        genome_all_posi = list(range(1,1+self.genome_len))
        modifible = list(set(genome_all_posi)-set(untouchble))
        r_positions = np.random.choice(modifible,size=2,replace=False)   
        return sorted(r_positions)

    def inversion(self, position1= None, position2=None): 
        '''
            input:
                inversion between position1 and position2
            rebust: 1. gene_len is a constant(l=1000)
        '''
        # Choose 2 position, position1<postiton2
        if position1 == None or position2 == None:
            r_position1,r_position2 = self.select_modify_position_inversion() #!!! not rubust need to change
            #print "delete position: "+str(r_position1)   
        else:
            # input position1/2 is for testing the code
            r_position1 = position1
            r_position2 = position2
        
        # inversion between r_posi1 and r_posi2
        Number_gene=0
        start_gene=[]
        end_gene=[]
        name_gene=[]
        orientation=[]
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
        # --- end inversion ----
        
        gene_to_modify = list(zip(name_gene, start_gene, end_gene,orientation))
        gene_to_modify.sort(key= lambda g:g[0])
        
        for i,g in enumerate([gn for gn in self.gene_list if gn.id in name_gene]):
            g.start = gene_to_modify[i][1]
            g.end = gene_to_modify[i][2]
            g.orientation = gene_to_modify[i][3]
        		
        return name_gene, start_gene, end_gene , orientation
        

'''
We can not draw genome in 5Bim's computer

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
def draw_genome(genome):
    
    g_list = genome.gene_list
    p_list = genome.prot_posi
    gdd = GenomeDiagram.Diagram('Test Diagram')
    gdt_features = gdd.new_track(1, greytrack=False) 
    gds_features = gdt_features.new_set()
    #Add three features to show the strand options,
    for g in g_list:
        feature = SeqFeature(FeatureLocation(g.start, g.end), strand=g.orientation)
        gds_features.add_feature(feature, name='gene'+str(g.id), label=True) #care for name
    
    for p in p_list:
        feature = SeqFeature(FeatureLocation(p, p+10))                                  
        gds_features.add_feature(feature,color='red',label=True,name='Protein position')
    
    # use feature _set   
    
    gdd.draw(format="circular", circular=True,circle_core=0.7, pagesize=(20*cm,20*cm),
             start=1, end=genome.genome_len) # careful for the length of genome
    gdd.write("GD_labels_default.pdf", "pdf")
'''


# check genome
'''
for g in gene_list:
    g.display()
'''

def check_choose_position(gn):
    '''choose postions_to_modify many times to see if they are inside the gene
       
        Arguments:
           gn(Genmoe): Genome instance
    '''
    g_list = gn.gene_list
    gene_range = [(g.start,g.end) for g in g_list]
    untouchble = []
    for r in gene_range:
        untouchble += range(r[0],r[1]+1)
    '''
    # test select_modify_position_inversion()
    test_position = []
    error = []
    for i in range(10000):        
        p = gn.select_modify_position_insert()
        #test_position.append(p)
        if  p in untouchble:
            error.append(p)
    # plt.hist(test_position)
    '''
   
    '''
    # test select_modify_position_inversion()
    test_position2 = []
    error2 = []
    for i in range(10000):        
        p = gn.select_modify_position_delete(60)
        #test_position2.append(p)
        if  p in untouchble:
            error2.append(p)
    return error
    # plt.hist(test_position2)
    '''
   
    test_position3 = []
    error3 = []
    for i in range(1000):        
        p = gn.select_modify_position_inversion()
        #test_position3.append(p)
        if  p[0] in untouchble or p[1] in untouchble:
            error3.append(p)
    return error3,test_position3
    
        
if __name__ == "__main__":    
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
            a.orientation = -1
        gene_list.append(a)    
    

        
    # ----- do some modification
    gn1 = Genome()
    gn1.gene_list=gene_list
    
    # ------ show genome------
    for g in gn1.gene_list:
        g.display()
    gn1.inversion(50, 16000)
    
    #gn1.insert(insert_len=6000)
    
    #gn1.delete(delete_len=3000)
    
    #gn1.gene_list.append(gene_list[1])
    
    for g in gn1.gene_list:
        g.display()
        
    # Bio-python draw genome, unfortunetely  we can't do it in 5bim
    #draw_genome(gn1)

  
