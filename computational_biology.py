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
            prot_posi(int list): list of protein position,with defaut value, attention
                        1.input manually might cause problem with diff prot. posi,
                        2.Also, assume that protine name = hns all the time
    '''
    def __init__(self):
        # each gene is a list [gene_id, start position, end positon, orientation]
        self.gene_list = []
        self.change_rate = 1.0
        self.generation = 0
        # defaut value is according to the TP 
        self.genome_len = 30000
        self.prot_posi = [1,3001,6001,9001,12001,15001,18001,21001,24001,27001] # initial, need to modify
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
        #print ("insert position: "+str(r_position))
        # find all genes after r_position, every gene + insert_len
        genes_to_modify = [g for g in self.gene_list if g.start > r_position]
        for g in genes_to_modify:
            g.start = g.start + insert_len
            g.end = g.end + insert_len
        
        #genome_len + insert_len
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
           #print ("delete position: "+str(r_position))
           
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
        
        # invert also protine list
        # ? nothing between posi1 and posi2, ? 1 prot between posi1 and posi2 --> OK
        self.prot_posi = self.modify_prot(method="inversion", posi1=r_position1,posi2=r_position2)
        
        # inversion between r_posi1 and r_posi2
        Number_gene=0
        start_gene=[]
        end_gene=[]
        name_gene=[]
        orientation=[]
        for g in self.gene_list:
            if (r_position1 < g.start and r_position2> g.end) : 
                Number_gene = Number_gene +1
                start_gene.append(g.start)
                end_gene.append(g.end)
                name_gene.append(g.id)
                orientation.append(g.orientation)
        
        # if Number_gene = 0, stop inversion process(nothing happens)
        if Number_gene == 0:
            return ("nothing changed")
        #print ('inversion between',r_position1,r_position2)         
        
        # sort gene list to inversion, by its start position
        g_list = list(zip(name_gene, start_gene, end_gene,orientation))        
        g_list.sort(key= lambda g:g[1])
        name_gene, start_gene, end_gene,orientation = [list(t) for t in zip(*g_list)]
        #print(name_gene, start_gene, end_gene,orientation)                
               
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
            self.inversion()
            mu = 3 #"inversion"
        
        return mu

    def loop_mutation(self, N) :  	
#N : number of mutations
	for i in range(N): 
		mut=self.mutation()
		if (mut ==1): 
			self.insert() 
		if (mut==2): 
			self.delete()
		if (mut==3): 
			self.inversion() 
		print mut		
		print self.gene_list 
		print self.prot_posi

	tousidentfile(self)
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
    
    # test select_modify_position_inversion()
    test_position = []
    error = []
    for i in range(10000):        
        p = gn.select_modify_position_insert()
        #test_position.append(p)
        if  p in untouchble:
            error.append(p)
    # plt.hist(test_position)
    
   
    
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
    
   
    test_position3 = []
    error3 = []
    for i in range(1000):        
        p = gn.select_modify_position_inversion()
        #test_position3.append(p)
        if  p[0] in untouchble or p[1] in untouchble:
            error3.append(p)
    return error3,test_position3
    
def check_mutation(gn,N = 1000):
    '''Test 3 methods of genome modification
    
        Arguments:
            gn(Genome): Genome instance
            N(int): Nbr of modifications on a genome(N generation)
    '''
    mutation_list = []
    for i in range(N):
        mu = gn.mutation()
        mutation_list.append(mu)
    return mutation_list

def check_modify_prot(gn, N=1000):
    '''Verify if Genome.modify_prot is correct
    
    '''
    #c = "inversion"
    
    '''
    # check modify_prot function
    for i in range(N):
        #c = np.random.choice(['insert','delete','inversion'])
        p1,p2 = 10,4000
        #print(p1,p2)
        gn.prot_posi = gn.modify_prot(c,posi1=p1,posi2=p2,length=60)
        print(gn.prot_posi)
    ''' 
    # check mutation with adding modify_prot
    print("original protine list:\n",gn.prot_posi)
    mu_list = []
    for i in range(N):
        mu = gn.mutation()
        #print(mu,gn.prot_posi)
        mu_list.append(mu)
    return mu_list
    
def tousidentfile(gn1):
	header=["##gff-version 3","#!gff-spec-version 1.20","#!processor NCBI annotwriter",
	"##sequence-region tousgenesidentiques 1 %d" % gn1.genome_len]
	f= open("./input/tousgenesidentiques1.gff","w+")
	for i in range(len(header)):
		f.write("%s\n" % header[i])
	f.write("tousgenesidentiques\tRefSeq\tregion\t1\t%d\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n" %  gn1.genome_len
	)
	for n in range(10):
		if (gn1.gene_list[n].orientation==1):
			gn1.gene_list[n].orientation='+'
		else :
			gn1.gene_list[n].orientation='-'
		base="tousgenesidentiques\tRefSeq\tgene\t%s\t%s\t.\t%s\t.\tID=g1;Name=g%s\n" % (gn1.gene_list[n].start,gn1.gene_list[n].end,gn1.gene_list[n].orientation,gn1.gene_list[n].id)
		f.write(base)
	f.close() 
	#for TSS.DAT
	header2=["TUindex","TUorient","TSS_pos","TSS_strength\n"]
	f2=open("./input/TSS.dat","w+")
	f2.write("%s" % header2[0])
	for i in range(1,len(header2)): 
		f2.write("\t%s" % header2[i])
	for n in range(10): 
		if (gn1.gene_list[n].orientation==1):
			gn1.gene_list[n].orientation='+'
            base="%s\t%s\t%s\t0.2\n" %(gn1.gene_list[n].id, gn1.gene_list[n].orientation, gn1.gene_list[n].start)
            f2.write(base)
        else :
			gn1.gene_list[n].orientation='-'
            base="%s\t%s\t%s\t0.2\n" %(gn1.gene_list[n].id, gn1.gene_list[n].orientation, gn1.gene_list[n].end)
            f2.write(base)	
	f2.close()
	#for TTS.dat
	header3=["TUindex","TUorient","TTS_pos","TTS_proba_off\n"]
	f3=open("./input/TTS.dat","w+")
	f3.write("%s" % header3[0])
	for i in range(1,len(header3)): 
		f3.write("\t%s" % header3[i])
	for n in range(10): 
		if (gn1.gene_list[n].orientation==1):
			gn1.gene_list[n].orientation='+'
            base="%s\t%s\t%s\t1.\n" %(gn1.gene_list[n].id, gn1.gene_list[n].orientation, gn1.gene_list[n].end)
            f3.write(base)  
		else :
			gn1.gene_list[n].orientation='-'
            base="%s\t%s\t%s\t1.\n" %(gn1.gene_list[n].id, gn1.gene_list[n].orientation, gn1.gene_list[n].start)
            f3.write(base)	
	f3.close()
	#for prot.dat 
	header4=["prot_name","prot_pos\n"]
	f4=open("./input/prot.dat","w+")	
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
    gn1.loop_mutation(10)
    #gn1.gene_list=gene_list
    #gn1.display_genome()
    
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

  
