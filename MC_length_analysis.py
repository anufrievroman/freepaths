import matplotlib.pyplot as plt
import numpy as np
import os

file_name_1='600nm.txt'
file_name_2='600nm_50K.txt'
file_name_3='600nm_100K.txt'

number_of_nodes=800
path_length_range=4000e-9
os.chdir('Length analysis')

def distribution_calculation(file_name):
    '''This function calculates the free path distribution from the raw free path data'''
    print '\n',file_name, 'is being processed'
    f = open(file_name,"r")
    data = np.loadtxt(f, dtype='float')
    f.close()
    distribution=np.zeros((number_of_nodes,2))
    distribution[:,0]=[i*path_length_range/(number_of_nodes) for i in range(number_of_nodes)]
    distribution[:,1]=[len(filter(lambda x: x>=j-0.5*path_length_range/number_of_nodes and x<j+0.5*path_length_range/number_of_nodes and x!=0, data)) for j in distribution[:,0]]
    return distribution

dist1=distribution_calculation(file_name_1)
dist2=distribution_calculation(file_name_2)
dist3=distribution_calculation(file_name_3)

plt.semilogy(dist1[:,0]*1e9,dist1[:,1],'b')
plt.semilogy(dist2[:,0]*1e9,dist2[:,1],'r')
plt.semilogy(dist3[:,0]*1e9,dist3[:,1],'g')
plt.ylabel('Number of phonons', fontsize=12)
plt.xlabel('Flight length along y (nm)', fontsize=12)
plt.savefig("Data analysis.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
plt.show()
