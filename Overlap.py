import numpy as np 
import cmath as cm 
import matplotlib.pyplot as plt 
from DDA import E_tot
#FIXME: Need to stop running whole DDA.py when importing E_tot -> see "if __name__ == "__main__"""

class target():
    #This will be the matrix that defines what field we want to see!
    def __init__(self, E_tot):
        self.E_tot = E_tot
        self.n_elements = len(E_tot)
        self.target = np.zeros(self.n_elements, dtype=np.complex_)

mytarg = target(E_tot)

for i in range(target(E_tot).n_elements):
    mytarg.target[i] =2+0.4*np.cos(i**2) + 0.4*1j*np.sin(i**2) # random target function

# now we have the matrix for the target function, we have to dot product with the field along the line
output = E_tot[:,400] 
print(output) # Output -> field along the line z = 400

normalisation_target = 0
normalisation_output = 0
#initialising normalisation constants for output and target vectors

for k in range(len(mytarg.target)):
    #Calculating normalisation constant for the vectors
    normalisation_target += np.real(mytarg.target[k])**2 + np.imag(mytarg.target[k])**2

normalisation_target = 1/np.sqrt(normalisation_target)

for l in range(len(E_tot)):
    normalisation_output += np.real(output[l])**2 + np.imag(output[l])**2

normalisation_output = 1/np.sqrt(normalisation_output)

target_hat = normalisation_target * mytarg.target
output_hat = normalisation_output * output

obj_value = np.abs(np.dot(output_hat, np.conj(target_hat)))
print(obj_value)

#TODO: make this all into a function so all we need to pass is where you want the output field to be measured at and the target function!


