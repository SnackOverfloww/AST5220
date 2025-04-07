
import numpy as np

n_k = 100
k_array = np.linspace(-4, 0, n_k)
print(k_array)

k_array_log = np.logspace(-4, 0, n_k)
print(k_array_log)

test = pow(10, k_array)
print(test)