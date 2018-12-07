# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:54:28 2018

@author: Lily
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
#sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle/') 
import simulation as sim 
from computational_biology import *
#import SamMayer.TCDS_v2.TCDS.plotting_v2 as plotting
#matplotlib.rcParams.update({'font.size': 13})

INI_file= 'params.ini'
output_dir='output'
output_dir_res = output_dir+"/all_res"


# ------ input Genome information -----------
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
# first generation
sim.start_transcribing(INI_file, output_dir)
#gn1.cal_fitness()


# mutation
fit_list = []
mutation_type = []
for i in range(20):
    gn_before=gn1 
    mut = gn1.mutation()
    #gn1.display_genome()
    mutation_type.append(mut)
    gn1.sort_by_start_posi()
    tousidentfile(gn1)
    try:
        print("mutation type =",mut)
        sim.start_transcribing(INI_file, output_dir)
    except:
        print ('mutation type', mut)
        gn1.display_genome()
        break
    f = gn1.call_fitness(gn1)
    f_before=gn1.call_fitness(gn_before)
    #monte carlo
    T=1/(i+1)
    r=gn1.monte_carlo(f_before,f,gn_before,T)
    print(r)
    if r==0: 
        gn1.gene_list=gn_before.gene_list
        print(r)
    
    fit_list.append(gn1.fitness)


np.savetxt('fitness.output', fit_list)    
np.savetxt('mutation.output', mutation_type)    

plt.plot(fit_list)
plt.show()  

