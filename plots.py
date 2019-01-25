#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 17:45:29 2018

@author: yyang1
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import re


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
    #plt.title('Initial parameters: T(0.1)_idL(60)_ratio(1)_polyNbr(6)')
    
    plt.show()


def plot_total_fit(fit_total, name_total,para,rep=5):
    f = np.array(fit_total)
    if para == 'ratio':
        n = list(set([n[23:-8] for n  in name_total]))
    elif para == 'idL':
         n = [n[13:-17] for n in name_total[:10]]
    elif para == 'T0':
        n = [n[:-27] for n in name_total[:7]]    
    elif para == 'poly':
        n = [n[34:-8] for n in name_total[:8]]
        
    colors = ['red','blue','green','black','yellow','pink','orange','gray','purple','brown']
    #n_l = name_total
    for j, name in enumerate(n):
        for i in range(rep):
            k = i*len(n) + j
            plt.plot(range(len(f[k])), f[k],c=colors[j])
    for j,name in enumerate(n):
        plt.plot([],[],label=name,c=colors[j])
    plt.xlabel('generations')
    plt.ylabel('fitness')
    plt.title('fitness evolution of different '+para)
    plt.legend()
    plt.show()
    #plt.savefig('./parameters_test/ratio/1')


def plot_init_para(fit_total, mutation_type, name_total, rep=5):
    
    for i in range(rep):
        plot_fitness(fit_total[i],mutation_type[i])

def time_reach_max(fit_total, name_total, rep = 5):
    tm = np.zeros([5,10])
    count = 0
    for i in range(5):
        for j in range(10):
            tm[i][j] = np.where(fit_total[count][:]==max(fit_total[count]))[0][0]
            count += 1
    return tm

def read_result(fit,mut,name,rep=5):
    nbr_para = int(len(name)/rep)
    fit_max = [[] for i in range(nbr_para)]
    
    count = 0
    for i in range(rep):
        for j in range(nbr_para):
            fit_max[j].append(max(fit[count]))
            count += 1
    #mu = np.zeros(len(mut),)
    return fit_max

def plot_max2(fit,mut,name,para,rep=5):
    f_max = read_result(fit,mut,name,rep)
    if para == 'T0':
        n = [n[:-27] for n in name[:7]] 
    elif para == 'idL':
        n = [n[13:-17] for n in name[:10]] 
    elif para == 'ratio':
        n = [n[23:-8] for n in name[:10]] 
    elif para == 'poly':
        n = [n[34:-8] for n in name[:8]] 
    
    plt.boxplot(f_max)
    plt.xticks(range(1,len(n)+1),n)
    plt.title('Test of different '+para)
    plt.xlabel(para)
    plt.ylabel('Max fitness value in 500 generations')
    plt.show()
    


#plot len
f_path1 = './result_data/len/'
fit = np.loadtxt(f_path1+'fit_total.output')
mut = np.loadtxt(f_path1+'./mutation_total.output')
name = np.genfromtxt(f_path1+'./name_total.output',dtype='str')
plot_max2(fit,mut,name,para='idL')
plot_total_fit(fit,name,para='idL',rep=5)

# ratio
f_path2 = './result_data/ratio/'
fit = np.loadtxt(f_path2+'fit_total.output')
mut = np.loadtxt(f_path2+'./mutation_total.output')
name = np.genfromtxt(f_path2+'./name_total.output',dtype='str')
plot_max2(fit,mut,name,para='ratio')
plot_total_fit(fit,name,para='ratio',rep=5)
tm = time_reach_max(fit,name)
plt.boxplot(tm)

# T0 test 2 
f_path4 = './result_data/T0/'
fit = np.loadtxt(f_path4+'fit_total.output')
mut = np.loadtxt(f_path4+'./mutation_total.output')
name = np.genfromtxt(f_path4+'./name_total.output',dtype='str')
plot_max2(fit,mut,name,para='T0')
plot_total_fit(fit,name,para='T0',rep=5)

# plot poly
# error fitness = -1
f_path5 = './result_data/poly/'
fit = np.loadtxt(f_path5+'fit_total3.output')
mut = np.loadtxt(f_path5+'./mutation_total3.output')
name = np.genfromtxt(f_path5+'./name_total3.output',dtype='str')
plot_max2(fit,mut,name,para='poly')

for f in fit :
    for i in range(len(f)):
        if f[i] == -1:
            f[i] = 0

plot_total_fit(fit,name,para='poly',rep=5)
#plot_total_fit(fit_total,name_total)
#plt.plot(fit_list)
#plt.show()  

# opt
f_path6 = './result_data/opt/'
fit = np.loadtxt(f_path6+'fit_total4.output')
mut = np.loadtxt(f_path6+'./mutation_total4.output')
name = np.genfromtxt(f_path6+'./name_total4.output',dtype='str')
plt.plot(fit[0])
plt.xlabel('generations')
plt.ylabel('fitness')
plt.title('fitness evolution with optimal parameters ')




