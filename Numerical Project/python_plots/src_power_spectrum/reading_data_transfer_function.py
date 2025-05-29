import numpy as np
import matplotlib.pyplot as plt

transfer_data = np.loadtxt("transfer_function.txt")

k_values = transfer_data[:,0]
data_k_1 = transfer_data[:,1] 
data_k_2 = transfer_data[:,2]
data_k_3 = transfer_data[:,3]
data_k_4 = transfer_data[:,4]
data_k_5 = transfer_data[:,5]



