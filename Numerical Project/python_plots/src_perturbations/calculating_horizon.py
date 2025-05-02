
import numpy as np
import matplotlib.pyplot as plt

#Imports the categorized data from the text file
exec(open("python_plots/src_perturbations/reading_data_perturbations.py").read())

# For the mode k = 0.001
eta_crossing_horizon_0_001 = ((1/0.001)*3.08567758e22)

i = 0
while eta_values_0_001_student[i] < eta_crossing_horizon_0_001:
    i+=1
bigger_than_horizon_0_001 = eta_values_0_001_student[i]
    
j = len(eta_values_0_001_student)
while eta_values_0_001_student[i] > eta_crossing_horizon_0_001:
    i-=1
smaller_than_horizon_0_001 = eta_values_0_001_student[i]

index_smaller_0_001 = (np.where(eta_values_0_001_student==smaller_than_horizon_0_001)[0][0]) 
index_bigger_0_001 = (np.where(eta_values_0_001_student==bigger_than_horizon_0_001)[0][0]) 
value_horizon_0_001 = (eta_values_0_001_student[index_smaller_0_001] + eta_values_0_001_student[index_bigger_0_001])/2 
x_value_horizon_0_001 = (x_values[index_smaller_0_001] + x_values[index_bigger_0_001])/2
x_value_uncertainty_horizon_0_001 = np.abs(x_values[index_smaller_0_001] - x_values[index_bigger_0_001])/2 

print("For the mode k = 0.001/Mpc:")
print("The value for n when kn = 1 is: ", str.format('{0:.2f}',value_horizon_0_001 / (60*60*24*365*(10**9)*2.99792458e8)), "Gyrs")
print("The value for x where kn = 1 is: ", 
    str.format('{0:.4f}', x_value_horizon_0_001), "+-", 
    str.format('{0:.4f}', x_value_uncertainty_horizon_0_001), "\n")


#========================================================================================================================

# For the mode k = 0.01
eta_crossing_horizon_0_01 = ((1/0.01)*3.08567758e22)

i = 0
while eta_values_0_01_student[i] < eta_crossing_horizon_0_01:
    i+=1
bigger_than_horizon_0_01 = eta_values_0_01_student[i]
    
j = len(eta_values_0_01_student)
while eta_values_0_01_student[i] > eta_crossing_horizon_0_01:
    i-=1
smaller_than_horizon_0_01 = eta_values_0_01_student[i]

index_smaller_0_01 = (np.where(eta_values_0_01_student==smaller_than_horizon_0_01)[0][0]) 
index_bigger_0_01 = (np.where(eta_values_0_01_student==bigger_than_horizon_0_01)[0][0]) 
value_horizon_0_01 = (eta_values_0_01_student[index_smaller_0_01] + eta_values_0_01_student[index_bigger_0_01])/2 
x_value_horizon_0_01 = (x_values[index_smaller_0_01] + x_values[index_bigger_0_01])/2
x_value_uncertainty_horizon_0_01 = np.abs(x_values[index_smaller_0_01] - x_values[index_bigger_0_01])/2 

print("For the mode k = 0.01/Mpc:")
print("The value for n when kn = 1 is: ", str.format('{0:.2f}',value_horizon_0_01 / (60*60*24*365*(10**6)*2.99792458e8)), "Myrs")
print("The value for x where kn = 1 is: ", 
    str.format('{0:.4f}', x_value_horizon_0_01), "+-", 
    str.format('{0:.4f}', x_value_uncertainty_horizon_0_01), "\n")

#===================================================================================================================================

# For the mode k = 0.1
eta_crossing_horizon_0_1 = ((1/0.1)*3.08567758e22)

i = 0
while eta_values_0_1_student[i] < eta_crossing_horizon_0_1:
    i+=1
bigger_than_horizon_0_1 = eta_values_0_1_student[i]
    
j = len(eta_values_0_1_student)
while eta_values_0_1_student[i] > eta_crossing_horizon_0_1:
    i-=1
smaller_than_horizon_0_1 = eta_values_0_1_student[i]

index_smaller_0_1 = (np.where(eta_values_0_1_student==smaller_than_horizon_0_1)[0][0]) 
index_bigger_0_1 = (np.where(eta_values_0_1_student==bigger_than_horizon_0_1)[0][0]) 
value_horizon_0_1 = (eta_values_0_1_student[index_smaller_0_1] + eta_values_0_1_student[index_bigger_0_1])/2 
x_value_horizon_0_1 = (x_values[index_smaller_0_1] + x_values[index_bigger_0_1])/2
x_value_uncertainty_horizon_0_1 = np.abs(x_values[index_smaller_0_1] - x_values[index_bigger_0_1])/2 

print("For the mode k = 0.1/Mpc:")
print("The value for n when kn = 1 is: ", str.format('{0:.2f}',value_horizon_0_1 / (60*60*24*365*(10**6)*2.99792458e8)), "Myrs")
print("The value for x where kn = 1 is: ", 
    str.format('{0:.4f}', x_value_horizon_0_1), "+-", 
    str.format('{0:.4f}', x_value_uncertainty_horizon_0_1), "\n")


#===================================================================================================================================

# For the mode k = 0.3

eta_crossing_horizon_0_3 = ((1/0.3)*3.08567758e22)

i = 0
while eta_values_0_3_student[i] < eta_crossing_horizon_0_3:
    i+=1
bigger_than_horizon_0_3 = eta_values_0_3_student[i]
    
j = len(eta_values_0_3_student)
while eta_values_0_3_student[i] > eta_crossing_horizon_0_3:
    i-=1
smaller_than_horizon_0_3 = eta_values_0_3_student[i]

index_smaller_0_3 = (np.where(eta_values_0_3_student==smaller_than_horizon_0_3)[0][0]) 
index_bigger_0_3 = (np.where(eta_values_0_3_student==bigger_than_horizon_0_3)[0][0]) 
value_horizon_0_3 = (eta_values_0_3_student[index_smaller_0_3] + eta_values_0_3_student[index_bigger_0_3])/2 
x_value_horizon_0_3 = (x_values[index_smaller_0_3] + x_values[index_bigger_0_3])/2
x_value_uncertainty_horizon_0_3 = np.abs(x_values[index_smaller_0_3] - x_values[index_bigger_0_3])/2 

print("For the mode k = 0.3/Mpc:")
print("The value for n when kn = 1 is: ", str.format('{0:.2f}',value_horizon_0_3 / (60*60*24*365*(10**6)*2.99792458e8)), "Myrs")
print("The value for x where kn = 1 is: ", 
    str.format('{0:.4f}', x_value_horizon_0_3), "+-", 
    str.format('{0:.4f}', x_value_uncertainty_horizon_0_3), "\n")

#===================================================================================================================================
