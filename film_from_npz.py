import sys
import numpy as np
import simulation as sim
import plotting as plotting
import matplotlib
matplotlib.rcParams.update({'font.size': 13})


sys.path.append("C:/Users/Lily/Documents/GitHub/Biologie-computationnelle")
INI_file="TCDS/params.ini"
output_dir="output"

output_dir_res = output_dir+"/all_res"

sigma_info = np.load(output_dir_res+"/save_sigma_info.npz")
RNAPs_info = np.load(output_dir_res+"/save_RNAPs_info.npz")

Barr_sigma_info = sigma_info["Barr_sigma_info"]
Dom_size_info = sigma_info["Dom_size_info"]

for i, Barr_sigma_val in enumerate(sigma_info["Barr_sigma_info"]):
	one_sigma_info = np.repeat(Barr_sigma_val, sigma_info["Dom_size_info"][i])
	RNAPs_pos_info = RNAPs_info["RNAPs_info"][:, 1, i]
	plotting.plot_mean_sigma_genes_v2(INI_file, one_sigma_info, RNAPs_pos_info)
