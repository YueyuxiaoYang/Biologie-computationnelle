# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:28:49 2018
@author: Yang Pauline Lilia
"""

import numpy as np
#import matplotlib.pyplot as plt

''' if you wan to import this file 
import sys
sys.path
# add the path of this file in sys so as to import it
sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle') 
from computational_biology import *
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
        attention:
            always start < end
    
    '''
    def __init__(self):
         # define defaut value
        self.id = 0
        self.start = -1 # start position should be between 1 and 30000
        self.end = 0 # end position should be [1:30000]
        self.orientation = 0 # '+':1, '-':0
        self.length = 999 # length of the gene is 1000 by defaut,and it is a constant 
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
            prot_posi(int list): list of protein position,with defaut value, attention
                        1.input manually might cause problem with diff prot. posi,
                        2.Also, assume that protine name = hns all the time
    '''

    def __init__(self):
        # each gene is a list [gene_id, start position, end positon, orientation]
        self.gene_list = []
        self.change_rate = 1
        self.generation = 0
        # defaut value is according to the TP 
        self.genome_len = 30000
        self.prot_posi = [1,3001,6001,9001,12001,15001,18001,21001,24001,27001] # initial, need to modify
        self.fitness = 0
        self.fit_exp = np.array([0.036,0.036,0.273,0.091,0.091,0.091,0.018,0.045,0.045,0.273])
        self.T0 = 0.00001 # how likely accept the 'worse' fitness value
    def display_genome(self):
        ''' display every gene in the genome 
        '''
        for g in self.gene_list:
            g.display()
            
    
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
        # find all genes after r_position, every gene + insert_len
        genes_to_modify = [g for g in self.gene_list if g.start > r_position]
        for g in genes_to_modify:
            g.start = g.start + insert_len
            g.end = g.end + insert_len
        
        self.genome_len += insert_len       
        # insert also protine list
        self.prot_posi = self.modify_prot(method="insert", posi1=r_position,length=insert_len)
        
    
    def delete(self, position=None, delete_len=60):
        '''delete a sequence in the genome
           Arguments:
               position (int): delete position, defaut random choose
               delete_len (int): the length of deletion, defaut=60
       '''
        if position == None:
           r_position = self.select_modify_position_delete(delete_len) 
        else:
           r_position = position
        # find all genes after r_position, every gene - delete_len
        genes_to_modify = [g for g in self.gene_list if g.start > r_position]   
        for g in genes_to_modify:
            g.start = g.start - delete_len
            g.end = g.end - delete_len
        # delete also protine list
        self.prot_posi = self.modify_prot(method="delete", posi1=r_position,length=delete_len)
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
        untouchble += self.prot_posi
        genome_all_posi = range(1,1+self.genome_len)
        modifible = list(set(genome_all_posi)-set(untouchble))
        r_position = np.random.choice(modifible)                
        return r_position
	
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
        # delete_posi - delete_len > prot_posi
        for p in self.prot_posi:
            untouchble += range(p,p+delete_len+1)
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
        untouchble += self.prot_posi
        genome_all_posi = list(range(1,1+self.genome_len))
        modifible = list(set(genome_all_posi)-set(untouchble))
        r_positions = np.random.choice(modifible,size=2,replace=False)   
        return sorted(r_positions)
    
    def modify_prot(self, method,posi1,posi2=None,length=None):
        ''' Modify the protine list 
            Arguments: 
                method(int): insert, delete or inversion
                posi1: insert/delete position, posi1 for inversion
                posi2: posi2 for inversion
            Rebust:
                1. every prot_posi produce same protine(ignore prot_id)
                2. input arguments not checked
        '''
        pl = self.prot_posi
        po = posi1
        if method == "insert":
            prot_list = [ps for ps in pl if ps < po]+[ps+length for ps in pl if ps > po]
        if method == "delete":
            prot_list = [ps for ps in pl if ps < po]+[ps-length for ps in pl if ps > po]
        if method == "inversion":
            prot_to_inv = [ps for ps in pl if ps > posi1 and ps < posi2]         
            pr_not_inv = [ps for ps in pl if ps not in prot_to_inv]
            #invert the prots
            prot_to_inv.append(posi2)
            pr_dist = np.array(prot_to_inv[::-1])
            pr_dist = pr_dist[:-1] - pr_dist[1:]
            p_after_inv = posi1 + np.add.accumulate(pr_dist)
            prot_list = list(p_after_inv) + pr_not_inv
        prot_list.sort()
        return prot_list
    
    def inversion_v2(self, posi1= None, posi2=None):
        if posi1 == None or posi2 == None:
            posi1,posi2 = self.select_modify_position_inversion() 
            #print "delete position: "+str(r_position1)   
        self.prot_posi = self.modify_prot(method="inversion", posi1=posi1,posi2=posi2)
        gene_l = [[g.id, g.start, g.end, g.orientation] for g in self.gene_list] # start, end and orientation
        gene_l.sort(key=lambda x:x[1]) # sort by start position
        g_to_inv = [g for g in gene_l if g[1] > posi1 and g[2] < posi2]         
        g_not_inv = [g for g in gene_l if g not in g_to_inv]
        id_l = [] # id list after sorting by start position
        s_e_l = [posi1,] # start and end positions
        ori_l = []
        for x in g_to_inv:
            id_l.append(x[0])
            s_e_l.append(x[1])
            s_e_l.append(x[2])
            ori_l.append(x[3])
        s_e_l.append(posi2)
        # inverse id;  orientation and  position will be update later
        g_to_inv = g_to_inv[::-1]
        #invert positions
        g_dist = np.array(s_e_l)
        g_dist = g_dist[1:] - g_dist[:-1]
        g_dist = g_dist[::-1]
        g_after_inv = posi1 + np.add.accumulate(g_dist) # inversion of positions complete
        g_after_inv = g_after_inv[:-1] # delete posi2
        # save start and end positions
        for i in range(len(g_to_inv)):
            g_to_inv[i][1] = g_after_inv[i*2]
            g_to_inv[i][2] = g_after_inv[i*2+1]
            g_to_inv[i][3] = g_to_inv[i][3] * (-1)
        gene_l_after_inv = g_to_inv + g_not_inv
        #print(gene_l_after_inv)
        # sort gene list aftere inversion by id
        gene_l_after_inv.sort(key=lambda x:x[0])
        #print(gene_l_after_inv)
        for i,g in enumerate([g for g in self.gene_list]):
            g.start = gene_l_after_inv[i][1]
            #print(i,g.start)
            g.end = gene_l_after_inv[i][2]
            g.orientation = gene_l_after_inv[i][3]
    
    def mutation(self):
        '''Randomly choose in/del or inversion for current genome, the result is next generation 
        '''
        r = self.change_rate
        if np.random.rand(1) < r/(r+1): # in/del
            if np.random.rand(1) < 0.5: # insert
                self.insert()
                mu = 1 #"insert"
            else: #delete
                self.delete()
                mu = 2 #"delete"
        else: # inversion
            self.inversion_v2()
            mu = 3 #"inversion"
        return mu
    
    def call_fitness(self,name=''):
        obs = np.genfromtxt('output'+name+'/save_tr_nbr.csv')
        obs = obs/sum(obs)
        fitness = np.exp(-1*sum(np.abs(obs-self.fit_exp)/self.fit_exp))
        self.fitness = fitness
        return self.fitness
    
    def sort_by_start_posi(self):
        gl = self.gene_list
        gl.sort(key = lambda x: x.start)    
    
    def monte_carlo(self,gn_new): 
        f_before = self.fitness
        f_new = gn_new.fitness
        if f_new<f_before : 
            #p= np.exp(-1/(1000*T))
            #k=1
            p= np.exp(-1*(f_before-f_new)/self.T0)##### mettre la bonne formule
            #p = 0.1
            r=np.random.binomial(1,p)
        else : 
            r=1
        return r 
                
            
            
        
    
'''
We can not draw genome in 5Bim's computer, however, you can try AWS-C9, a online virtual IDE
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

    
def tousidentfile(gn1,name=''):
	''' This function permits to write in files
	'''
    header=["##gff-version 3","#!gff-spec-version 1.20","#!processor NCBI annotwriter",
	"##sequence-region tousgenesidentiques 1 %d" %gn1.genome_len]
    f= open("./input"+name+"/gff.gff","w+")
    for i in range(len(header)):
        f.write("%s\n" % header[i])
    f.write("tousgenesidentiques\tRefSeq\tregion\t1\t%d\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n" %gn1.genome_len)
    for n in range(10):
        if (gn1.gene_list[n].orientation==1):
            base="tousgenesidentiques\tRefSeq\tgene\t%s\t%s\t.\t%s\t.\tID=g1;Name=g%s\n" %(gn1.gene_list[n].start,gn1.gene_list[n].end,"+",gn1.gene_list[n].id-1)            
        else :
            base="tousgenesidentiques\tRefSeq\tgene\t%s\t%s\t.\t%s\t.\tID=g1;Name=g%s\n" %(gn1.gene_list[n].end,gn1.gene_list[n].start,"-",gn1.gene_list[n].id-1)
        f.write(base)
    f.close() 
	#for TSS.DAT
    header2=["TUindex","TUorient","TSS_pos","TSS_strength\n"]
    f2=open("./input"+name+"/TSS.dat","w+")
    f2.write("%s" % header2[0])
    for i in range(1,len(header2)): 
        f2.write("\t%s" % header2[i])
    for n in range(10):
        if (gn1.gene_list[n].orientation==1):
            orient="+"
            base="%s\t%s\t%s\t1.\n" %(gn1.gene_list[n].id-1, orient, gn1.gene_list[n].start)
        else :
            orient="-"
            base="%s\t%s\t%s\t1.\n" %(gn1.gene_list[n].id-1, orient, gn1.gene_list[n].end)
        f2.write(base)	
    f2.close()
	#for TTS.dat
    header3=["TUindex","TUorient","TTS_pos","TTS_proba_off\n"]
    f3=open("./input"+name+"/TTS.dat","w+")
    f3.write("%s" % header3[0])
    for i in range(1,len(header3)): 
        f3.write("\t%s" % header3[i])
    for n in range(10): 
        if (gn1.gene_list[n].orientation==1):
            orient="+"
            base="%s\t%s\t%s\t1.\n" %(gn1.gene_list[n].id-1, orient, gn1.gene_list[n].end)
        else :
            orient="-"
            base="%s\t%s\t%s\t1.\n" %(gn1.gene_list[n].id-1, orient, gn1.gene_list[n].start)
        f3.write(base)	
    f3.close()
	#for prot.dat 
    header4=["prot_name","prot_pos\n"]
    f4=open("./input"+name+"/prot.dat","w+")	
    f4.write("%s" % header4[0])
    for i in range(1,len(header4)): 
        f4.write("\t%s" % header4[i])	
    for n in range(10): 
        base="hns\t%s\n" %(gn1.prot_posi[n])
        f4.write(base)
    f4.close()	
	
        
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
    tousidentfile(gn1)
    
    # verify mutation
    #cl = check_mutation(gn1)
    #plt.hist(cl)
   # mu_list = check_mutation(gn1)
   # np.hist(mu_list)
    # ------ show genome------

    # --- fichier 
    #gn1.insert()
    #gn1.gene_list.sort(key= lambda g:g.start) # sort the gene list by its start posi
    #tousidentfile(gn1)
        
    # Bio-python draw genome, unfortunetely  we can't do it in 5bim
    #draw_genome(gn1)
