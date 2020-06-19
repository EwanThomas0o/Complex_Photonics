import numpy as np
import matplotlib as mpl
import cmath as cm
import scipy as sp
import matplotlib.pyplot as plt
import os 


class dipole():
    def __init__(self, r, polarisability):                               
        self.r = r 
        self.polarisability = polarisability
        self.P = 0                         

    # Arguments:
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

	# Arguments:
	# 	modk is the modulus of the wavevector
	#	amplitude is the amplitude of the incoming electric plane wave
	#	t is the time that we freeze the system at 
	#	w is the angular frequency of the wave

	def E_at_given_R(self, r):
		E = cm.rect(self.Amp, (np.dot(self.k, r) - self.w * self.t))
		# We can use cm.rect() to act as the complex exponential. First argument is the prefactor to the exp(), second is what's exponentiated
		return E


def createddipoles(positions):
	#TODO: make is so that this function reads a CSV file and then puts the dipoles where they belong
	volume = np.zeros((500,500,500))
	# Inistialising total volume simulated
	dipoles = []
	# Initialising array to store dipole information


	for i in range(len(positions)):
		r = positions[i]
		dipoles.append(dipole(r, 10))
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
	# Inverse of the interaction matrix -> this is the slow step in DDA
	Polarisations = np.matmul(A_inverse, E_at_dipoles)
	# Obtaining exact polarisations of dipoles.

	
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
						E_wave += Ajk * dic["dipoles"][k].P
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

def detector(incoming_wave, y, z, dic):
	# Function designed to calculate the field profile "far" away from the dipole arrangement
	plot_vector = np.zeros(500, dtype = np.complex_)
	# Field values are stored in a 1D column vector as complex numbers 
	for i in range(0,500):
		point = [i,y,z]
		E_wave = incoming_wave.E_at_given_R(point)
		# Find the value of the field due to inout wave
		for k in range(len(dic["dipoles"])):
			rjk = np.linalg.norm(np.subtract(point, dic["dipoles"][k].r))
			# Modulus of distance from point to kth dipole

			if rjk != 0:
				# As long as we're not on the dipole (remember, E diverges here)
				Ajk = cm.rect( (1/rjk), (incoming_wave.modk * rjk))
				E_wave += Ajk * dic["dipoles"][k].P
				plot_vector[i] = E_wave
				# Adding the contribution due to the dipoles emitted radiation					
				# Giving plot area a value of E everywhere
	return plot_vector



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
		# Sum of real and imaginary parts squared -> defintion of the complex dot product
	normalisation_target = 1/np.sqrt(normalisation_target)
	# Normalisation constant for target vector

	output_hat = normalisation_output * output
	target_hat = normalisation_target * target
	# Normalising each vector

	obj_value = np.abs(np.dot(output_hat, np.conj(target_hat)))
	# How much one vector (and thus field) "overlaps" another
	return obj_value

# Some defined values ------------------------------------------
# /////////////////////////////////////////////////////////////

y = 200
z = 400
t = 15
incoming_wave = Wave(2*np.pi/10, 2, t)
positions_1 = [[100,200,200],[300,200,200]]
positions_2 = [[100,200,200],[293,200,213]]
dic_1 = createddipoles(positions_1)
plot_area_1 = dic_1["volume"][:, y, :].copy().astype(np.complex_)
polarisations_1 = calculateP(incoming_wave, dic_1)
dic_2 = createddipoles(positions_2)
polarisations_2 = calculateP(incoming_wave, dic_2)
# Setting some values so we can plot a picture!


#E_tot_1 = CalculateE(incoming_wave, plot_area_1, y, dic_1)

ff_1 = detector(incoming_wave, y, z, dic_1)
# Farfield of first arangement of dipoles
ff_2 = detector(incoming_wave, y, z,dic_2)
# Farfield of second arangement of dipoles

#print(detector(positions_1, incoming_wave, y, 400, dic_1))
print(overlap(ff_2,ff_1))
#print(dic_1["dipoles"][0].r, dic_1["dipoles"][1].r)

def newmovedipole(dic, ff_1, ff_2, incoming_wave):
	# This is a function uses gradient descent in order to find the optimal position of dipoles in order to match a target field
	# This is so much faster that the gradient descent script! Takes seconds rather than around 10 minutes!
	for j in range(len(dic["dipoles"])):
		primary_dipole = dic["dipoles"][j]

		x0 = primary_dipole.r[0]
		z0 = primary_dipole.r[2]
		# Storing the initial positions of the jth dipoles
		ovr = [overlap(ff_1, ff_2)]  
		# Will store the overlaps of the field due to the change in position of the dipoles 
		pos = []
		# Will store positions of the dipole when moved in a circle 
		pos.append(primary_dipole.r)
		# This represents original position's overlap and postion
		
		t = 0
		while t <= 2*np.pi:
			primary_dipole.r = [x0 + 1*np.cos(t), 200, z0 + 1*np.sin(t)] 
			# Moving in a cicle around each dipole
			point = primary_dipole.r
			pos.append(point)
			
			
			# Need to recalculate polarisation of dipole in this position before we calc the overlap P = alphaE
			
			E_wave = incoming_wave.E_at_given_R(point)
			for k in range(len(dic["dipoles"])):
				if j != k:
					rjk = np.linalg.norm(np.subtract(point, dic["dipoles"][k].r))
					Ajk = cm.rect( (1/rjk), (incoming_wave.modk * rjk))
					E_wave += Ajk * dic["dipoles"][k].P	
					# 0th order approximation employed here. We assume that the other dipoles don't really change in repsonse to the jth dipole moving
				else:
					pass
			# Calculating the Electric field at everypoint on the circle to find the polarisation of the dipole.
			# If we know the polarisations we can find the Electric field anywhere! Including on the "detector"
			primary_dipole.P = primary_dipole.polarisability * E_wave
			
			new_field = detector(incoming_wave, 200, 400, dic)
			# Calculating the field due to the new postion
			ovr.append(overlap(ff_2, new_field))
			# Then find the overlap 

			t += 0.1
			# Itterate around the circle 

		index = ovr.index(max(ovr))
		# Finding the position of the maximum overlap
		primary_dipole.r = pos[index]
		# Changing the postion of the jth dipole that corresponds to the greatest value of overlap
		calculateP(incoming_wave, dic) #Resolve DDA to obtain exact Polarisations
		print(primary_dipole.r, ovr[index]) 
		
	
	return newmovedipole(dic, detector(incoming_wave, y, z, dic), ff_2, incoming_wave  ) 
	# Call itteratively 

# TODO: migrate this function to gradient_descent_rehaul.py

print(newmovedipole(dic_1, ff_1, ff_2, incoming_wave)["dipoles"][1].r)



# Some plotting commands -----------------------------------
# //////////////////////////////////////////////////////////

# vals = np.abs(ff_1)   #What we've got here is the instantaeous electric field, not the actual field 
# vals_2d = np.real(E_tot_1)  #E_tot is the vector that needs to be dot producted with the TARGET vector.
# x = np.linspace(0,499,num=500)
# z = np.linspace(0,499, num=500)


# fig = plt.figure(1)
# plt.plot(x,vals)
# plt.xlabel('Position along line')
# plt.ylabel('Magnitude of Electric field')
# plt.show()


# fig = plt.figure(2)
# ax = plt.imshow(vals_2d, cmap = 'cool', interpolation = 'nearest')
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

def createanimation(ani_vals_dip_1, ani_vals_dip_2):
    fig, ax = plt.subplots()
    
    def animate(i):
        ax.clear()
        ax.set(xlim=(100, 310), ylim=(190, 210))
        ax.plot(ani_vals_dip_1[i][0], ani_vals_dip_1[i][2],'o')
        ax.plot(ani_vals_dip_2[i][0], ani_vals_dip_2[i][2],'o')
        ax.plot(100, 200, 'o')
        plt.annotate("Stationary Dipole", (100,200))
        ax.plot(293,213, '+')
        plt.annotate("Target Location", (293,213))
        plt.xlabel("Z Axis")
        plt.ylabel("X Axis")
        

        

    animation = ani.FuncAnimation(fig, animate, 10, interval = 20)

    plt.show()
    return animation

#///////////////////////////////////////////////////////////////////////////


