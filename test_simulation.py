# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:54:28 2018

@author: Lily
"""

sys.path.append('/home/yyang1/Bureau/Biologie-computationnelle/SamMayer/TCDS-v2-master/TCDS/') 
import simulation as sim 


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
