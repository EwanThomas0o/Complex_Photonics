from DDA import overlap, farfield, dipole, Wave
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# This script will make a copy of the positions 2 vector in the DDA folder and then it evaluate the overlap integral
#   of two functions when dipole is moved in 4 different directions, then it will alter the copied vector to move
#   in the direction where the overlap value is the largest, it will do this util it reaches a minimum 

positions_1 = [[100,200,200],[301,200,203]]
positions_2 = [[100,200,200],[297,200,198]]  #This is the dipole we will move

y = 200
t = 15
incoming_wave = Wave(2*np.pi/10, 2, t)


target = farfield(positions_1, y ,t ,incoming_wave) 
out = farfield(positions_2, y ,t ,incoming_wave)        
# These functions yield what the field looks like at z = 400


original_overlap = overlap(target,out)
# This is the overlap between the dipoles initial arangement and the target arangement

ani_vals = []
# Simple storage of the location of the dipoles to be animated.

ani_overlap = []
# imple storage of the values of the overlap to be displayed in the title.

def movedipoles(original_overlap, positions_2, target, y, t, incoming_wave, ani_vals, ani_overlap): 
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
        ovr.append(overlap(target, farfield((new_positions[i]), y, t, incoming_wave)))
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
    ani_vals.append(np.asarray(info["positions"][no][1]))
    ani_overlap.append(mx)

    if mx == original_overlap:
        return ani_vals , ani_overlap 
        # If the max overlap is the same as the overlap from the previous step then we've hit a minimum! (No guarantee it's global)
    
    else: 
        return movedipoles(mx, info["positions"][no], target, y, t, incoming_wave, ani_vals, ani_overlap )
        # If the max overlap is different then we continue moving in the direction of increasing overlap!



print(ani_vals)
print(ani_overlap)

def createanimation(ani_vals, ani_overlap):
    #Animates the dipoles moving to an optimal position
    fig, ax = plt.subplots()
    
    def animate(i):
        ax.clear()
        # Clears plot after every iteration
        ax.set(xlim=(100, 310), ylim=(190, 210))
        ax.plot(ani_vals[i][0], ani_vals[i][2],'o')
        # Plot the position of the dipole that moves, this is the one that changes (see the i)
        plt.annotate("Moving Dipole", (ani_vals[i][0], ani_vals[i][2]))
        # Label that moves with the dipole
        ax.plot(100, 200, 'o')
        # Plot of the stationary dipole
        plt.annotate("Stationary Dipole", (100,200))
        # Label for the stationary dipole
        ax.plot(301,203, '+')
        # Plot of the target position
        plt.annotate("Target Location", (301,203))
        # Label for the target postion
        plt.xlabel("Z Axis")
        plt.ylabel("X Axis")
        plt.title('Overlap at step = {}'.format(ani_overlap[i]))
        # Puts the value of the overlap of the two functions in the title
        

        

    animation = ani.FuncAnimation(fig, animate, 10, interval = 20)

    plt.show()
    return animation
    
ani_vals, ani_overlap = movedipoles(original_overlap, positions_2, target, y, t, incoming_wave, ani_vals, ani_overlap)
animation =createanimation(ani_vals, ani_overlap)

os.system('say "your program has finished"')
# Audio Notification





