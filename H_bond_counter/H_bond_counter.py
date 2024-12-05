import sys
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np

#*XDATCAR.xlsx file as input to the script
input_file = sys.argv[1]
dataset = pd.read_excel(input_file)
dataset = dataset.fillna(0)

#lattice parameters in angstrom
a = (dataset.iloc[[0], :]).values.tolist()[0][0]
b = (dataset.iloc[[1], :]).values.tolist()[0][1]
c = (dataset.iloc[[2], :]).values.tolist()[0][2]

#atoms and number of atoms
atom = (dataset.iloc[[3], :]).values.tolist()[0]
number_of_atom = (dataset.iloc[[4], :]).values.tolist()[0]
atom_list = sum([[s] * n for s, n in zip(atom, number_of_atom)], [])

#coordinates
X = np.array((dataset.iloc[5:,[0]]).values.tolist())
Y = np.array((dataset.iloc[5:,[1]]).values.tolist())
Z = np.array((dataset.iloc[5:,[2]]).values.tolist())

#frames
number_of_frames = math.floor(len(X)/sum(number_of_atom))

if (len(X) != len(Y)):
    print("Check Input File. Coordinates missing/extra!")

if (len(X) != len(Z)):
    print("Check Input File. Coordinates missing/extra!")

#raising errors for missing frame data
if ((len(X)/sum(number_of_atom)) > number_of_frames): 
    print("Missing data in the last frame!")
    sys.exit()

#reference atom number for tracking the h bonds
index_reference_atom = int(sys.argv[2]) - 1

if (index_reference_atom > sum(number_of_atom) - 1):
    print("Reference Atom Index Out Of Range")
    sys.exit()

#check if the reference atom can form h bond or not
if atom_list[index_reference_atom] not in ['F', 'O', 'N']:
    print("The reference atom can not form hydrogen bonds")
    sys.exit()

#specify the cut-off distance for h bonds in Angstroms
h_bond_cutoff = 2.0
covalent_bond_cutoff = 1.2

# Identify all H can form h-bond
h_bonded = np.zeros(len(atom_list))

for hydrogens in range(sum(number_of_atom)):
    if atom_list[hydrogens] == 'H':
        x_atom1, y_atom1, z_atom1 = X[hydrogens][0], Y[hydrogens][0], Z[hydrogens][0]
        bond_distance = np.array([])
        distance = []
        for all_atoms in range(sum(number_of_atom)):
            periodic_distance = [[[0 for k in range(3)] for j in range(3)] for i in range(3)]
            x_atom2, y_atom2, z_atom2 = X[all_atoms][0], Y[all_atoms][0], Z[all_atoms][0]
            for m in range(-1, 2):
                for n in range(-1, 2):
                    for o in range(-1, 2):
                        periodic_distance[m+1][n+1][o+1] = math.dist((x_atom1*a, y_atom1*b, z_atom1*c), ((x_atom2 + m)*a, (y_atom2 + n)*b, (z_atom2 + o)*c))
            distance += [min(min(min(periodic_distance)))]
            if distance[all_atoms] == 0:
                distance[all_atoms] = 100
        bond_distance = np.append(bond_distance,[min(distance)])
        if atom_list[distance.index(min(distance))] in ['F', 'O', 'N']:
            h_bonded[hydrogens] = 1

#initializing distances and counters
periodic_distance = [[[0 for k in range(3)] for j in range(3)] for i in range(3)]
h_bond_counter = np.zeros(number_of_frames)
avg_h_bond_counter = np.zeros(number_of_frames)

#calculate distances and use cut-off distances to count H-bondings
for frames in range(number_of_frames):
    x_reference_atom, y_reference_atom, z_reference_atom = X[frames*sum(number_of_atom) + index_reference_atom][0], Y[frames*sum(number_of_atom) + index_reference_atom][0], Z[frames*sum(number_of_atom) + index_reference_atom][0]
    distance = [0]*(sum(number_of_atom))
    #print(h_bond_counter[frames])
    for index_target_atom in range(sum(number_of_atom)):
        x_target_atom, y_target_atom, z_target_atom = X[frames*sum(number_of_atom) + index_target_atom][0], Y[frames*sum(number_of_atom) + index_target_atom][0], Z[frames*sum(number_of_atom) + index_target_atom][0]
        for m in range(-1, 2):
            for n in range(-1, 2):  
                for o in range(-1, 2):
                    periodic_distance[m+1][n+1][o+1] = math.dist((x_reference_atom*a, y_reference_atom*b, z_reference_atom*c), ((x_target_atom + m)*a, (y_target_atom + n)*b, (z_target_atom + o)*c))
        distance[index_target_atom] = min(min(min(periodic_distance)))
        if distance[index_target_atom] == 0:
            distance[index_target_atom] = 100
        if (covalent_bond_cutoff <= distance[index_target_atom] <= h_bond_cutoff) and (h_bonded[index_target_atom] == 1):
             h_bond_counter[frames] += 1  
    avg_h_bond_counter[frames] = sum(h_bond_counter)/(frames+1)

#calculate average and standard deviation of hydrogen bonds
std = np.std(h_bond_counter)
avg = avg_h_bond_counter[-1]
stats = np.array([[avg, std]])

#plotting and printing output
fig = plt.figure(figsize=(10, 6))
plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.style'] = 'normal'
plt.plot(h_bond_counter, color = "black")
plt.plot(avg_h_bond_counter, "red", label = "Average", linewidth = 3)
plt.xlabel('Number of Frames')
plt.xlim([0, number_of_frames])
plt.ylabel('Number of hydrogen bonds \n for the reference atom')
plt.legend(loc="upper right")
output_file = sys.argv[3]
plt.savefig(output_file+'.png')


#print and save output
print(f'Avg Number of H Bonding: {np.round(avg,4)} \nStandard deviation: {np.round(std,4)}')
np.savetxt(f"{output_file}.txt", stats, header="Average Standard Deviation", comments='', fmt='%1.4e')
