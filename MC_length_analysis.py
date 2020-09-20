import matplotlib.pyplot as plt
import numpy as np
import time

path = ["Length analysis/Different segments/400nm.txt",
        "Length analysis/Different segments/600nm.txt",
        "Length analysis/Different segments/800nm.txt"]

number_of_nodes = 400
path_length_range = 10000e-9


def distribution_calculation(filename, data_range, number_of_nodes):
    '''This function calculates distribution of numbers in a given file'''
    data = np.loadtxt(filename)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:,0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:,1], _ = np.histogram(data, number_of_nodes, range=(0, data_range))
    return distribution


for filename in path:
    dist = distribution_calculation(filename, path_length_range, number_of_nodes)
    plt.plot(dist[:,0]*1e9, dist[:,1], 'b')


plt.ylabel('Number of phonons', fontsize=12)
plt.xlabel('Flight length along (nm)', fontsize=12)
plt.savefig("Data analysis.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
plt.show()

print ("Mean free path is", np.mean(data)*1e9, "nm")
