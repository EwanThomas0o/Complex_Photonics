import numpy as np
import matplotlib as mpl
import cmath as cm
import matplotlib.pyplot as plt
import scipy as sc


class dipole():
    def __init__(self, r, polarisability):                               
        self.r = r 
        self.polarisability = polarisability
        self.P = 0                         

    # Variables of Class:
    #   r is dipole position
    #   polarisability is used to calculate the matrix A
    #   P is the Polarisation of each dipole used to calculate reemitted field

class Wave():
	def __init__(self, modk, amplitude, time):
		self.modk = modk
		self.k = modk * np.asarray([0,0,1])
		self.Amp = amplitude
		self.t = time
		self.w = self.modk * 1
		# c = 1 here (can use SciPy if physical values wanted)

	# Variables of Class:
	# 	modk is the modulus of the wavevector
	#	amplitude is the amplitude of the incoming electric plane wave
	#	t is the time that we freeze the system at 
	#	w is the angular frequency of the wave

	def E_at_given_R(self, r):
		E = cm.rect(self.Amp, (np.dot(self.k, r) - self.w * self.t))
		# We can use cm.rect() to act as the complex exponential. First argument is the prefactor to the exp(), second is what's exponentiated
		return E

	# Arguements:
	# 	r is the position of evaluation

def createddipoles(positions):
	#TODO: make is so that this function reads a CSV file and then puts the dipoles where they belong
	volume = np.zeros((500,500,500))
	# Initialising total volume simulated
	dipoles = []
	# Initialising array to store dipole information


	for i in range(len(positions)):
		r = positions[i]
		dipoles.append(dipole(r, 50))
		volume[r[0]][r[1]][r[2]] = 1
	
	dic = {"dipoles" : dipoles, "volume": volume}
	# Dictionary that stores the dipoles and total volume. Location of dipoles are denoted by a 1 in the volume matrix
	# dipoles stores information on the dipoles, like their position and polarisability
	return dic

def calculateA(incoming_wave, dic):
	#Finding the A matrix to then find the polarisations of each dipole
	n = len(dic["dipoles"])
	# Number of dipoles
	A = np.empty((n,n), dtype = np.complex_)
	# Initiallising an empty matrix with complex data type.

	for j in range(n):
		primary_dipole = dic["dipoles"][j]
		# Focus on the of the dipoles then consider interaction with all others, then switch dipole and do the same.
		for k in range(n):
			if j != k:
				secondary_dipole = dic["dipoles"][k]
				
				rjk = np.linalg.norm( np.subtract(primary_dipole.r, secondary_dipole.r) ) 
				A[j][k] = cm.rect( (1 / rjk) , (incoming_wave.modk * rjk) )
				#See paper by Draine 1994 for equans
			else: 
				A[j][j] = 1 / (primary_dipole.polarisability)

	return A

def calculateP(incoming_wave, dic):
	A = calculateA(incoming_wave, dic)
	#Finding the interaction matrix for the dipoles
	n = len(dic["dipoles"])
	E_at_dipoles = np.zeros(n, dtype = np.complex_)
	# Initialising array for the value of E at the dipoles
	
	for i in range(n):
		E_at_dipoles[i] = incoming_wave.E_at_given_R(dic["dipoles"][i].r)
		#Contains information about phase and amplitude at the point of the dipoles.
	
	A_inverse = np.linalg.inv(A)
	# Slow step as A is rather large

	Polarisations = np.matmul(A_inverse, E_at_dipoles)

	
	for i in range(n):
		dic["dipoles"][i].P = Polarisations[i]
		# Assigning the polarisation values to each of the dipoles
	
	return Polarisations
	# Polarisation values are important as they allow us to calculate the field due to the dipole at any point 

def CalculateE(incoming_wave, plot_area, y, dic):
	for i in range(plot_area.shape[0]): 
		# Plot area is 500x500 so this func is going through each point and is summing electric field from all contributions at that point
		for j in range(plot_area.shape[1]):

			if plot_area[i][j] == 1:
				# We don't evaluate as r is zero and this electric field diverges
				pass
			
			else:
				point = [i , y, j]
				E_wave = incoming_wave.E_at_given_R(point)
				# Finds electric field at every point in the grid due to the incoming wave only
				

				for k in range(len(dic["dipoles"])):
					rjk = np.linalg.norm(np.subtract(point, dic["dipoles"][k].r))
					# Modulus of distance from point to kth dipole

					if rjk != 0:
						# As long as we're not on the dipole (remember, E diverges here)
						Ajk = cm.rect( (1/rjk), (incoming_wave.modk * rjk))
						E_wave -= Ajk * dic["dipoles"][k].P
						# Adding the contribution due to the dipoles emitted radiation
						plot_area[i][j] = E_wave
						# Giving plot area a value of E everywhere



	for dipole in dic["dipoles"]:
			r = dipole.r
			for i in range(-5, 5):   #Diameter is 10
				for j in range(-5, 5):
					plot_area[r[0] + i][r[2] + j] = 0
	return plot_area
	# Done to minimse the huge values of E near the dipole
	# TODO: introduce cross sectional scattering in order obtain actual values of E on the dipoles

def farfield(positions, y ,t ,incoming_wave):
	# A function simply designed to retrieve the far field of a given arrangement of dipoles
	dic = createddipoles(positions)
	polarisations = calculateP(incoming_wave,dic)
	plot_area = dic["volume"][:, y, :].copy().astype(np.complex_)

	E_tot = CalculateE(incoming_wave, plot_area, y, dic)

	return E_tot[:,400] 
	# Very slow function! See DDA_rehaul.py for better farfield function (named detector in DDA_rehaul.py)
	

def overlap(target, output):
	# This fucntion finds "how similar" the far field of two arangements of dipoles are.
	# Target and output are vectors that give the values of electric field along a line.
	# These vectors are obtained from the function farfield()
	normalisation_target = 0
	normalisation_output = 0
	#initialising normalisations constants for vectors output and target
	for i in range(len(output)):
		normalisation_output += np.real(output[i])**2 + np.imag(output[i])**2

	normalisation_output = 1/np.sqrt(normalisation_output)
	# Normalisation constant for output vector

	for i in range(len(target)):
		normalisation_target += np.real(target[i])**2 + np.imag(target[i])**2

	normalisation_target = 1/np.sqrt(normalisation_target)
	# Normalisation constant for target vector

	output_hat = normalisation_output * output
	target_hat = normalisation_target * target
	# Normalising each vector

	obj_value = np.abs(np.dot(output_hat, np.conj(target_hat)))
	# How much one vector overlaps another
	return obj_value


y = 200
t = 15
incoming_wave = Wave(2*np.pi/10, 2, t)
positions_1 = [[100,200,200],[300,200,200]]
positions_2 = [[100,200,200],[290,200,200]]
dic_1 = createddipoles(positions_1)
plot_area_1 = dic_1["volume"][:, y, :].copy().astype(np.complex_)
polarisations_1 = calculateP(incoming_wave, dic_1)
# Setting some values to plot a picture!


E_tot_1 = CalculateE(incoming_wave, plot_area_1, y, dic_1)

ff_1 = farfield(positions_1, y ,t ,incoming_wave)
ff_2 = farfield(positions_2, y ,t ,incoming_wave)

print(overlap(ff_1,ff_2))

# Some plotting commands -----------------------------------
# //////////////////////////////////////////////////////////



# vals = np.real(E_tot_1)   #What we've got here is the instantaeous electric field, not the actual field 
# vals_2d = np.abs(farfield(positions_1, y, t, incoming_wave))  #E_tot is the vector that needs to be dot producted with the TARGET vector.
# x = np.linspace(0,499,num=500)
# z = np.linspace(0,499, num=500)


# fig = plt.figure(1)
# plt.plot(x,vals_2d)
# plt.xlabel('Position along line')
# plt.ylabel('Magnitude of Electric field')
# plt.show()


# fig = plt.figure(2)
# ax = plt.imshow(vals, cmap = 'cool', interpolation = 'nearest')
# plt.xlabel('Z-axis')
# plt.ylabel('X-axis')
# plt.title('Instantaneous Power of Electric Field at t = 15s')
# plt.colorbar(ax)
# plt.show()

#vals_2 = np.real(E_tot_2)   #What we've got here is the instantaeous electric field, not the actual field 
#vals_2d_2 = np.abs(E_tot_2[:,400])  #E_tot is the vector that needs to be dot producted with the TARGET vector.


# fig = plt.figure(3)
# plt.plot(x,vals_2d_2)
# plt.xlabel('Position along line')
# plt.ylabel('Magnitude of Electric field')
# plt.show()


# fig = plt.figure(4)
# ax = plt.imshow(vals_2, cmap = 'cool', interpolation = 'nearest')
# plt.xlabel('Z-axis')
# plt.ylabel('X-axis')
# plt.title('Instantaneous Power of Electric Field at t = 15s')
# plt.colorbar(ax)
# plt.show()


#///////////////////////////////////////////////////////////////////////



	




