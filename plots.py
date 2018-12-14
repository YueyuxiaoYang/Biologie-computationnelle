#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 17:45:29 2018

@author: yyang1
"""
import numpy as np
import matplotlib as plt
import sys

sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle/')
sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle/parameters_test/ratio/')

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

# read results
fit = np.loadtxt('./parameters_test/ratio/fit_total.output')


#plot_total_fit(fit_total,name_total)
#plt.plot(fit_list)
#plt.show()  
