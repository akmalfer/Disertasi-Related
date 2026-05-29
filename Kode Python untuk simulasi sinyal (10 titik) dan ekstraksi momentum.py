import os
import shutil
import numpy as np

# Get the directory where your Python script is located
python_dir = os.path.dirname(os.path.abspath(__file__))

# Define the correct path for the CSV file
momentum_file = os.path.join(python_dir, "10p_signal1jtestdatawmvis_falsedecays.csv")

# Open the file in write mode and add a header
with open(momentum_file, "w") as storage:
    storage.write("MSc, px_muons, py_muons, pz_muons, px_electrons, py_electrons, pz_electrons, px_neutrinos, py_neutrinos, mvis\n")

msc = np.linspace(201, 449, 10)

for m in msc:
    runback = "run_var_msc__testdata_falsedecays" + str(m) + "_v3"

    # Move to MadGraph directory and generate events
    os.chdir("/home/akmal/madgraph/MG5_aMC_v3_5_7/testscalartamudecay/")
    os.system('cp Cards/param_card_default.dat Cards/param_card.dat')
    
    lines = open("./Cards/param_card.dat", 'r').readlines()
    lines[31] = "   35 " + str(m) + " # MSc\n"
    
    with open("./Cards/param_card.dat", 'w') as out:
        out.writelines(lines)
    
    '''
    os.system(f'./bin/madevent generate_events {runback}')
    os.system(f'../Delphes/root2lhco Events/{runback}/tag_1_delphes_events.root Events/{runback}/tag_1_delphes_events.lhco')
    os.system(f'rm ./Events/{runback}/tag_1_delphes_events.root')
    '''
    
    # Move to MadAnalysis directory
    os.chdir("/home/akmal/madanalysis/madanalysis5-main/emunu_pdata/Build/")

    # Create input file for MadAnalysis
    with open("list.txt", "w") as my_file:
        my_file.write(f"/home/akmal/madgraph/MG5_aMC_v3_5_7/testscalartamudecay/Events/{runback}/tag_1_delphes_events.lhco")

    # Run MadAnalysis
    os.system('./MadAnalysis5job list.txt')

     # Read and process momentum data
    if os.path.exists("output_momentum.txt"):
        with open("output_momentum.txt", "r") as output_data:
            for line in output_data:
                values = line.strip().split()  # Split line into list of values
                if len(values) == 9:
                    values = list(map(float, values))  # Convert all values to floats
    
                    # Unpack the 9 values correctly
                    px_muons, py_muons, pz_muons, px_electrons, py_electrons, pz_electrons, px_neutrinos, py_neutrinos, vis_mass= values

                    with open(momentum_file, "a") as storage:
                        storage.write(f"{m}, {px_muons}, {py_muons}, {pz_muons}, {px_electrons}, {py_electrons}, {pz_electrons}, {px_neutrinos}, {py_neutrinos}, {vis_mass}\n")


        print(f"Momentum data stored for MSc = {m}")

print("All mass points processed successfully!")
