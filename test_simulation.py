# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:54:28 2018

@author: Lily
"""

import sys
#import matplotlib.pyplot as plt
import numpy as np
#sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle/') 
import simulation as sim 
from computational_biology import *
import copy
#import SamMayer.TCDS_v2.TCDS.plotting_v2 as plotting
#matplotlib.rcParams.update({'font.size': 13})



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
'''
def plot_fitness(fit_list, mutation_type):
    f = np.array(fit_list)
    #markers = {1:'v', 2:'P', 3:'o'}
    #color = {1:'red', 2:'blue', 3:'green'}

    insert = [i for i,t in enumerate(mutation_type )if t==1]
    delete = [i for i,t in enumerate(mutation_type )if t==2]
    inv = [i for i,t in enumerate(mutation_type )if t==3]
    
    plt.plot(range(len(fit_list)),fit_list)    
    plt.scatter(insert,f[insert],marker='v',c='red',label='ins')
    plt.scatter(delete,f[delete],marker='P',c='blue',label='del')
    plt.scatter(inv,f[inv],marker='o',c='green',label='inv')
    plt.legend()
    plt.show()


def plot_total_fit(fit_total, name_total):
    f = np.array(fit_total)
    n = name_total
    for i,name in enumerate(n):
        plt.plot(range(len(f[i])), f[i],label=name)
    plt.legend()
    plt.savefig('1')
'''
# ----- do some modification

#the prob. ratio between insertion/deletion and inversion.
#id_inv_ratio = [1/100,1/50, 1/10, 1, 10/1, 50/1, 100/1, 500/1, 1000/1, 10000/1] #
idL_list = [20, 40, 60, 80, 100, 120, 140, 200, 400, 800]

ratio = 1
T0 = 0.000001
#idL = 60
name = str(1) # paralell
iter_num = 500 # number of generations for one genome
rep_num = 5
fit_total = []
name_total = []
mutation_total = []

for rep in range(rep_num):
    for idL in idL_list:
        gn1 = Genome()
        gn1.T0 = T0
        gn1.idL = idL
        gn1.gene_list=gene_list
        gn1.change_rate = ratio
        
        INI_file= 'params'+name+'.ini'
        output_dir='output'+name
        output_dir_res = output_dir+"/all_res"
        
        tousidentfile(gn1,name)
        # first generation
        sim.start_transcribing(INI_file, output_dir)
        #gn1.cal_fitness()
        
        
        # mutation
        fit_list = []
        mutation_type = []
        g_list = []
        for i in range(iter_num):
            gn_new = copy.deepcopy(gn1)
            mut = gn_new.mutation()
            #gn1.display_genome()
            mutation_type.append(mut)
            gn_new.sort_by_start_posi()
            tousidentfile(gn_new,name)
            try:
                print("mutation type =",mut)
                sim.start_transcribing(INI_file, output_dir)
            except:
                print ('mutation type', mut)
                gn_new.display_genome()
                break
            f = gn_new.call_fitness(name)
            #f_before=gn.call_fitness(gn_before)
            #monte carlo
            T=1/(i+1)
            r=gn1.monte_carlo(gn_new)
            print("iter: idL=%s, ratio=%s, rep=%s, generation=%s"%(idL,ratio,rep,i))
            if r==1: # changer 
                gn1 = gn_new
                print('gnome changed')
            
            fit_list.append(gn1.fitness)
        fit_total.append(fit_list)
        mutation_total.append(mutation_type)
        name_total.append("T(%s)_idL(%s)_ratio(%s)_rep(%s)" %(T0,idL,ratio,rep))
            #plot_fitness(fit_list,mutation_type)
            #g_list.append(gn1.gene_list)
        #save files
        #np.savetxt("./parameters_test/ratio/fitness_T(%s)_idL(%s)_ratio(%s)_rep(%s).output" %(T0,idL,ratio,rep), fit_list)    
        #np.savetxt("./parameters_test/ratio/mutation_T(%s)_idL(%s)_ratio(%s)_rep(%s).output" %(T0,idL,ratio,rep), mutation_type)    
        np.savetxt("./parameters_test/len/fit_total.output", fit_total) 
        np.savetxt("./parameters_test/len/mutation_total.output", mutation_total)
        np.savetxt("./parameters_test/len/name_total.output", name_total, delimiter=" ", fmt="%s")
        #np.savetxt("./parameters_test/ratio/genelist_T(%s)_idL(%s)_ratio(%s)_rep(%s).output" %(T0,idL,ratio,rep), g_list)    
        

