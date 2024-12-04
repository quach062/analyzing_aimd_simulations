import sys
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np

#*XDATCAR.xlsx file as input to the script
input_file = sys.argv[1]
dataset = pd.read_excel(input_file)
dataset = dataset.fillna(0)

a = (dataset.iloc[[0], :]).values.tolist()[0][0]
b = (dataset.iloc[[1], :]).values.tolist()[0][1]
c = (dataset.iloc[[2], :]).values.tolist()[0][2]

atom = (dataset.iloc[[3], :]).values.tolist()[0]
number_of_atom = (dataset.iloc[[4], :]).values.tolist()[0]
X = np.array((dataset.iloc[5:,[0]]).values.tolist())
Y = np.array((dataset.iloc[5:,[1]]).values.tolist())
Z = np.array((dataset.iloc[5:,[2]]).values.tolist())

#raising errors for invalid input file
if (len(X) != len(Y)):
    print("Check Input File. Coordinates missing/extra!")
    sys.exit()

if (len(X) != len(Z)):
    print("Check Input File. Coordinates missing/extra!")
    sys.exit()

#frames in dynamic simualtion
number_of_frames = math.floor(len(X)/sum(number_of_atom))

#raising errors for missing frame data
if ((len(X)/sum(number_of_atom)) > number_of_frames):  
    print("Missing data in the last frame!") 
    sys.exit()

ref_frame = int(sys.argv[2])  #input frame number starting from 0 - num_frame-1
if ref_frame >= number_of_frames:
    print("Reference frame index is out of range!")
    sys.exit()

periodic_distance = [[[0 for k in range(3)] for j in range(3)] for i in range(3)]
total_distance = avg_total_distance = np.array([])

# RMSD with a specified ref frame
total_ref_distance = np.array([])
avg_total_ref_distance = np.array([])
for frames in range(number_of_frames - 1):
    #print(frames)
    ref_frame_distance = 0
    for atom in range(sum(number_of_atom)):
        distance = 0
        x_atom1, y_atom1, z_atom1 = X[(ref_frame)*sum(number_of_atom) + atom][0], Y[(ref_frame)*sum(number_of_atom) + atom][0], Z[(ref_frame)*sum(number_of_atom) + atom][0]
        x_atom2, y_atom2, z_atom2 = X[(frames + 1)*sum(number_of_atom) + atom][0], Y[(frames + 1)*sum(number_of_atom) + atom][0], Z[(frames + 1)*sum(number_of_atom) + atom][0]
        for m in range(-1, 2):
            for n in range(-1, 2):
                for o in range(-1, 2):
                    periodic_distance[m+1][n+1][o+1] = math.dist((x_atom1*a, y_atom1*b, z_atom1*c), ((x_atom2 + m)*a, (y_atom2 + n)*b, (z_atom2 + o)*c))
        distance = min(min(min(periodic_distance))) 
        ref_frame_distance += distance**2
    #print(ref_frame_distance)
    total_ref_distance = np.append(total_ref_distance, [math.sqrt((ref_frame_distance)/sum(number_of_atom))])
    avg_total_ref_distance = np.append(avg_total_ref_distance, [math.sqrt(sum(np.square(total_ref_distance))/len(total_ref_distance))])

fig = plt.figure(figsize=(10, 4))
plt.plot(total_ref_distance, color = "black")
plt.plot(avg_total_ref_distance, "red", label = "Average", linewidth = 3)
plt.xlabel('Number of Frames')
plt.xlim([0, number_of_frames])
plt.ylabel(r'Root Mean Square Displacement ($\AA$)' + f'\nwith respect to reframe number {ref_frame}')
plt.legend(loc="upper right")
output_file = sys.argv[3]
plt.savefig(output_file+'.png') # Adjust format of figure here

avg = avg_total_ref_distance[-1]
stats = np.array([[avg]])

print(f'Avg Distance: {np.round(avg,4)}')

# Save results in text file
np.savetxt(f"{output_file}.txt", stats, header="Average", comments='', fmt='%1.4e')