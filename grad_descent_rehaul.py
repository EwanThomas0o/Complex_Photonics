from DDA_rehaul import *
import numpy as np
import os

# This script will make a copy of the positions 2 vector in the DDA folder and then it evaluate the overlap integral
#   of two functions when dipole is moved in 4 different directions, then it will alter the copied vector to move
#   in the direction where the overlap value is the largest, it will do this util it reaches a minimum 

positions_1 = [[100,200,200],[301,200,203]]
positions_2 = [[100,200,200],[297,200,198]]  #This is the dipole we will move
dic_1 = createddipoles(positions_1)
polarisations_1 = calculateP(incoming_wave, dic_1)
dic_2 = createddipoles(positions_2)
polarisations_2 = calculateP(incoming_wave, dic_2)
y = 200
z = 400
t = 15
incoming_wave = Wave(2*np.pi/10, 2, t)


target = detector(positions_1,incoming_wave, y, z, dic_1) 
out = detector(positions_2,incoming_wave, y, z, dic_1)        
# These functions yield what the field looks like at z = 400


original_overlap = overlap(target,out)
#This is the overlap between the dipoles initial arangement and the target arangement


def movedipoles(original_overlap, positions_2, target, y, t, incoming_wave): 
    # This is gonna be a slow function as it will need to call farfield() from DDA four times, with two dipoles that'll take 100s
    #   per check of the local area...
    #   All in all takes around 400s for example provided below

    #List of positions made by moving the dipole, including original
    ovr = []
    #List of values of overlaps formed by the new fields created
    new_positions = []

    for i in range(0,5):
        new_positions.append(positions_2)
    new_positions = np.asarray(new_positions)

    step_size = 1

    new_positions[0][1][0] = new_positions[0][1][0] + step_size
    new_positions[1][1][2] = new_positions[0][1][2] + step_size
    new_positions[2][1][0] = new_positions[2][1][0] - step_size
    new_positions[3][1][2] = new_positions[3][1][2] - step_size

    # We now have 4 slightly altered versions of positons_1 along with the original, this is done the explore the values of overlap 
    #   in the vicinity of the moving dipole, we move in the direction of increasing overlap (as this means fields are more alike)

    for i in range(len(new_positions)-1):
        # Don't want to recalculate the overlap for the last pair of positions as this is given by "original_overlap" 
        ovr.append(overlap(target, detector((new_positions[i]), incoming_wave, y, z, dic)) 
        #TODO: WE NEED TO EVALUATE THE FIELD
        # AROUND THE DIPOLES TO FIND THE POLARISATIONS AGAIN, THIS IS WHAT DAVE MEANT WHEN HE WAS TALKING ABOUT FINDING THE FIELD IN 
        # THE REGION AROUND THE DIPOLES!!!
        # Calculating the new overlaps for each position
    ovr.append(original_overlap)
    # In the first iteration, this appends the original overlap. In every itteration after it appends the max overlap from the previous step 

    info = {"positions":new_positions, "overlaps":ovr}
    # dictionary that hold list of positions and also a list of the overlaps. Indexes correspond to respective values.
    # i.e. info["positions"][i] has overlap info["overlaps"][i] with the target function
    mx = max(info["overlaps"])
    # max value of all overlaps
    no = info["overlaps"].index(mx)
    # index of mx tells us which direction to move in

    print(info["positions"][no][1], mx)
    # Display postion of moving dipole and overlap for each step to terminal

    if mx == original_overlap:
        return 
        # If the max overlap is the same as the overlap from the previous step then we've hit a minimum! (No guarantee it's global)
    
    else: 
        return movedipoles(mx, info["positions"][no], target, y, t, incoming_wave )
        # If the max overlap is different then we continue moving in the direction of increasing overlap!

movedipoles(original_overlap, positions_2, target, y, t, incoming_wave)


os.system('say "your program has finished"')
# Notification




