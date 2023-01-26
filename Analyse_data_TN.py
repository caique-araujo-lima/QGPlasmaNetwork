import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous
from scipy.optimize import curve_fit
import math
import warnings
import sys
import os
import os.path
dir = os.getcwd()

N = 1000
d = 3
aA = 3.00
aG = 1.00
w0 = 1.00
eta = 2
a = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
a_interaction = 1.0

plt.style.use('dark_background')
cores = ['red', 'blue', 'yellow', 'green', 'violet', 'white']

bins = 100


def analise(N, d, aA, aG, w0, eta, a, a_interaction, bins):

    # Folder Path
    path1 = r"C:\Users\caiqu\OneDrive\Área de Trabalho\Dados"
    path2 = '\\N_%d_d_%d_aA_%.2f_aG_%.2f_w0_%.2f_eta_%.2f_a_%.6f_a_interaction_%.3f' %(N, d, aA, aG, w0, eta, a, a_interaction)
    path = path1 + path2
    energy_list = np.array([])
    print(path)
    # Change the directory
    os.chdir(path)
    # iterate through all file
    for file in os.listdir():
        # Check whether file is in text format or not
        if file.endswith(".txt"):
            file_path = f"{path}\{file}"
            print(file_path)
            # call read text file function
            tabela = np.loadtxt(file_path)
            energy_list = np.append(tabela[:,0], energy_list)

    energy_list = np.array( [ num for num in energy_list if num >= 0 ] )
    start, finish=min(energy_list), max(energy_list)

    bins_list=list(np.logspace(np.log10(start), np.log10(finish), num=bins))
    hist, edges=np.histogram(energy_list, bins=bins_list, density=True)

    p_e= hist/np.sum(hist)
    bins_midlist=edges[:-1]


    file = 'N_%d_d_%d_aA_%.2f_aG_%.2f_w0_%.2f_eta_%.2f_a_%.3f_a_interaction_%.3f_histogram_bins_%d.txt' %(N, d, aA, aG, w0, eta, a, a_interaction, bins)
    name = os.path.join(path, file)

    f = open(name, "w")
    for i in range (N):
        f.write(str(bins_midlist[i]) + '\t' + p_e(label_list[i]) + '\n')
    f.close()
    return 'DONE'

bins_midlist[i], p_e[i] = analise(N, d, aA, aG, w0, eta, a[0], a_interaction, bins)
