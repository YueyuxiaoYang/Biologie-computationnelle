# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:54:28 2018

@author: Lily
"""

import sys
import matplotlib
import numpy as np
#sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle/') 
import simulation as sim 
#import SamMayer.TCDS_v2.TCDS.plotting_v2 as plotting
#matplotlib.rcParams.update({'font.size': 13})

INI_file= 'params.ini'
output_dir='output'

output_dir_res = output_dir+"/all_res"
sim.start_transcribing(INI_file, output_dir)

'''
try:
    #output_dir=sys.argv[2]
    
except:
    print("sim goes wrong")
'''
'''
sigma_info = np.load(output_dir_res+"/save_sigma_info.npz")
RNAPs_info = np.load(output_dir_res+"/save_RNAPs_info.npz")

Barr_sigma_info = sigma_info["Barr_sigma_info"]
Dom_size_info = sigma_info["Dom_size_info"]

for i, Barr_sigma_val in enumerate(sigma_info["Barr_sigma_info"]):
	one_sigma_info = np.repeat(Barr_sigma_val, sigma_info["Dom_size_info"][i])
	RNAPs_pos_info = RNAPs_info["RNAPs_info"][:, 1, i]
	plotting.plot_mean_sigma_genes_v2(INI_file, one_sigma_info, RNAPs_pos_info)
'''


'''
if __name__ == '__main__':
    try:
        INI_file="params.ini" #sys.argv[1]           # e.g "../analysis_scripts/example/params.ini"
        # First simulation
        sim.start_transcribing(INI_file, "output")
        # or you can specify the output path
        #start_transcribing(INI_file, output_path)

        # Resuming
        # uncomment this two lines if you want to resume the simulation
        # 'first_output_path' means the path from which the script will get the npz files
        #first_output_path=sys.argv[2]        # e.g "../analysis_scripts/example/first_output"
        #start_transcribing(INI_file, first_output_path, resume=True)
    except configparser.NoSectionError:
        print("Error ! Please check the path to the paraneters files !")
        sys.exit(1)
    except PermissionError:
        print("Permission denied ! Please check the directory to the output files !")
        sys.exit(1)
    except (FileNotFoundError, NameError):
        print("Error ! Please check the directory to the output files !")
        sys.exit(1)
'''
