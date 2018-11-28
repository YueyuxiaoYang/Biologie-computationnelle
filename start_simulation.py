import sys
sys.path.append('C:/Users/Lily/Documents/GitHub/Biologie-comptationnelle')
import simulation as sim
INI_file="params.ini"
try:
    output_dir="output"
    sim.start_transcribing(INI_file, output_dir)
except:
    sim.start_transcribing(INI_file)
