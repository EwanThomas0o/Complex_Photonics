import numpy as np 

original_positions = [[100,200,200],[300,200,200]]

new_positions = []

for i in range(0,5):
    new_positions.append(original_positions)


new_positions = np.asarray(new_positions)


new_positions[0][1][0] = new_positions[0][1][0] + 5
new_positions[1][1][2] = new_positions[0][1][2] + 5
new_positions[2][1][0] = new_positions[2][1][0] - 5
new_positions[3][1][2] = new_positions[3][1][2] - 5

print(new_positions)

info = 






