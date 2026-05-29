import shutil
import numpy as np
import os


output_file = "testdataWW0j1mruwmvisv1wjets.txt"

n = 100  # Number of iterations

with open(output_file, "w") as results:  # Open output file for writing
    for m in range(1, n + 1):  # Loop from 1 to n
        runback = f"BackWWhtamutestdata_r{m}_v6_1mrun"

        # Move to MadGraph directory and generate events
        os.chdir("/home/akmal/madgraph/MG5_aMC_v3_5_7/pp_ww_emu_CMS/")
        '''
        os.system(f'./bin/madevent generate_events {runback}')
        os.system(f'../Delphes/root2lhco Events/{runback}/tag_1_delphes_events.root Events/{runback}/tag_1_delphes_events.lhco')
        os.system(f'rm ./Events/{runback}/tag_1_delphes_events.root')
        '''
        # Move to MadAnalysis directory
        os.chdir("/home/akmal/madanalysis/madanalysis5-main/emunu_pdata_withjets/Build/")

        # Create input file for MadAnalysis
        with open("list.txt", "w") as my_file:
            my_file.write(f"/home/akmal/madgraph/MG5_aMC_v3_5_7/pp_ww_emu_CMS/Events/{runback}/tag_1_delphes_events.lhco")

        # Run MadAnalysis
        os.system('./MadAnalysis5job list.txt')

        # Read the output file
        with open('output_momentum.txt', "r") as file:
            content = file.readlines()

        # Read and process momentum data
            with open("output_momentum.txt", "r") as output_data:
                for line in output_data:
                    values = line.strip().split()  # Split line into list of values
                    if len(values) == 12:
                        values = list(map(float, values))  # Convert all values to floats
    
                        # Unpack the 8 values correctly
                        px_muons, py_muons, pz_muons, px_electrons, py_electrons, pz_electrons, px_neutrinos, py_neutrinos, px_jets, py_jets, pz_jets, vis_mass= values

                        with open(output_file, "a") as storage:
                            storage.write(f"{px_muons}, {py_muons}, {pz_muons}, {px_electrons}, {py_electrons}, {pz_electrons}, {px_neutrinos}, {py_neutrinos}, {px_jets}, {py_jets}, {pz_jets}, {vis_mass}\n")