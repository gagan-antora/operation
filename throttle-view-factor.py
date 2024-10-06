import numpy as np
import matplotlib.pyplot as plt

# Define the geometry

receiver_length = 10  # Length of the receiver
receiver_position = (0.0,receiver_length/2)

panel_length = 2
panel_position = (4,receiver_length/2)  # X-position of the blocking panel (center of the panel)
panel_angle = 75 #Between 0.1 and 89.9 where 90 is fully open
panel_angle_range = np.arange(-89.9, 89.9, 2)



#diffuser_position = 10.0  # Position of the diffuse surface (emitter)
diffuse_surface_length = 8  # Length of the diffuse surface (emitter)
diffuser_position = (panel_position[0]+ diffuse_surface_length/2,receiver_length/2 - panel_length/2)  # Position of the receiver

diffuser_start_x = diffuser_position[0] - diffuse_surface_length/2
diffuser_start_y = diffuser_position[1]
diffuser_end_x = diffuser_position[0] + diffuse_surface_length/2
diffuser_end_y = diffuser_position[1]

barrier_length = diffuse_surface_length
barrier_position = (diffuser_position[0], diffuser_position[1] + panel_length)

# Simulation parameters
num_particles = 100000
counter = 0
list= []

def generate_old__particle():
    position = np.random.uniform(0, diffuse_surface_length)  # Random position on the emitte
    rand_cos = np.random.uniform(0, 1) # Distribution has to follow lambers law for diffuse surfaces
    angle = np.arccos(rand_cos)-np.pi/2
    return position, angle

def generate_particle_path(midpoint, length):
    x_position = np.random.uniform(midpoint[0]-length/2, midpoint[0]+length/2)  # Random position on the emitte
    rand_cos = np.random.uniform(-1, 1) # Distribution has to follow lambers law for diffuse surfaces
    angle = np.arccos(rand_cos)
    #print(angle)
    slope = np.tan(angle)  # Slope of the trajectory (dy/dx)
    y_inter = midpoint[1] - slope *x_position
    return slope, y_inter

# Check if the particle hits a vertical panel
def hits_vertical_panel(particle_slope, particle_intercept):
    global counter
    y_at_panel = particle_slope * (panel_position[0]) + particle_intercept  # y-position when x = panel_position
    if (panel_position[1] + panel_length/2) > (y_at_panel) > (panel_position[1] - panel_length/2) :
         # Particle hits the panel
        return True
    return False
def hits_horizontal_panel(particle_slope, particle_intercept):
    global counter
    x_at_barrier = (barrier_position[1] - particle_intercept ) /particle_slope
    if (diffuser_start_x) <= (x_at_barrier) <= (diffuser_end_x) :
        counter += 1
        return True
    return False

def hits_slanted_panel(particle_slope, particle_intercept, panel_angle):
    global counter
    panel_angle = np.radians(90-panel_angle)
    panel_start_x = panel_position[0] - (panel_length/2) * np.cos(panel_angle)
    panel_start_y = panel_position[1] - (panel_length/2) * np.sin(panel_angle)
    panel_end_x = panel_position[0] + (panel_length/2) * np.cos(panel_angle)
    panel_end_y = panel_position[1] + (panel_length/2) * np.sin(panel_angle)

    
    panel_slope = np.tan(panel_angle)
    panel_intercept = panel_start_y - panel_slope * panel_start_x
     
    # Solve for x_intersection:
    if particle_slope != panel_slope:  # Ensure the slopes aren't equal (no intersection)
        
        x_intersection = (panel_intercept - particle_intercept) / (particle_slope - panel_slope)
        y_intersection = particle_slope * (x_intersection) + particle_intercept
        y_intersection2 = panel_slope * (x_intersection) + panel_intercept
        if np.abs(y_intersection - y_intersection2) < 0.0001:
            pass
        # Check if the intersection point is within the bounds of the slanted panel
        if panel_start_x <= x_intersection <= panel_end_x or panel_start_y <= y_intersection <= panel_end_y:
            
            return True  # The particle hits the panel
        
    return False  # No hit

def hits_receiver(particle_slope, particle_intercept):

    y_at_receiver = particle_slope * (receiver_position[0]) + particle_intercept  # y-position at receiver x = 1.0
    if 0 <= y_at_receiver <= receiver_length:
        return particle_slope, y_at_receiver  # The particle hits the receiver
    return False

def view_factor_parallel_plates(w_i, w_j, L):
    # Calculate the normalized widths
    W_i = w_i / L
    W_j = w_j / L
    F_ij = (np.sqrt((W_i + W_j)**2 + 4) - np.sqrt((W_i - W_j)**2 + 4)) / (2 * W_i)
    return F_ij
def view_factor_perpendicular_plates(emmitor, panel):
    F_ij = (1 + (panel/emmitor) - ( 1 + (panel/emmitor)**2 )**0.5)/2
    return F_ij

slanted_panel_vf_range = []
slanted_reciever_vf_range = []

for panel_angle in panel_angle_range:
    # Simulation loop
    vertical_reciever_hits = 0
    vertical_panel_hits = 0
    slanted_reciever_hits = 0
    slanted_panel_hits = 0


    for _ in range(num_particles):
        particle_slope, particle_intercept = generate_particle_path(diffuser_position, diffuse_surface_length)  # Generate random particle emission
        if not hits_horizontal_panel(particle_slope, particle_intercept):
            if particle_slope < 0: #Don't include ones that are pointing away
            
                # Check if it hits horizontal the panel
                if hits_vertical_panel(particle_slope, particle_intercept):  # Check if it hits vertical the panel
                    vertical_panel_hits += 1
                else:
                    if hits_receiver(particle_slope, particle_intercept):  # Check if it hits the receiver
                        vertical_reciever_hits += 1
                
                if hits_slanted_panel(particle_slope, particle_intercept, panel_angle):  # Check if it hits the slanted panel
                    slanted_panel_hits += 1
                else:
                    hit = hits_receiver(particle_slope, particle_intercept)
                    if hit is not False:  # Check if it hits the receiver
                        list.append(hit)  # Store the particle path for later analysis
                        
                        slanted_reciever_hits += 1

    # Estimate view factor
    vertical_panel_vf = vertical_panel_hits / num_particles
    vertical_reciever_vf = vertical_reciever_hits / num_particles
    slanted_panel_vf = slanted_panel_hits / num_particles
    slanted_reciever_vf = slanted_reciever_hits / num_particles
    slanted_panel_vf_range.append(slanted_panel_vf)
    slanted_reciever_vf_range.append(slanted_reciever_vf)

analytical_vertical_panel_vf = view_factor_perpendicular_plates(diffuse_surface_length, panel_length)
analytical_vertical_reciever_vf = analytical_vertical_panel_vf * view_factor_parallel_plates(diffuse_surface_length, receiver_length, np.abs(receiver_position[0]-panel_position[0]))

print(counter)
print("--------------------------------------------------------")
print("Estimated blocks-vertical throttle VF:", slanted_panel_vf)
print("Estimated blocks-HX Box VF:", slanted_reciever_vf)
print("--------------------------------------------------------\n")
print("Max blocks-throttle VF:", analytical_vertical_panel_vf)
print("Max blocks-HX Box VF:", analytical_vertical_reciever_vf)

plt.plot(panel_angle_range, slanted_panel_vf_range, label='blocks-throttle', marker='o' )
plt.plot(panel_angle_range, slanted_reciever_vf_range, label='blocks-HX', marker='o')
plt.legend()
plt.show()
# print('panel start: ', panel_start_x, panel_start_y)
# print('panel end: ', panel_end_x, panel_end_y)
# print('diffuser start: ', diffuser_start_x, diffuser_start_y)
# print('diffuser end: ', diffuser_end_x, diffuser_end_y)

# print(view_factor_parallel_plates(8,8,2))

def plot_light_trajectory(slope_intercept_list):
    # Define x values between 0 and 10
    x = np.linspace(0, 30, 100)
    
    # Create a plot
    plt.figure(figsize=(8, 6))
    
    # Loop over each slope and intercept tuple in the list
    for slope, intercept in slope_intercept_list:
        # Calculate y values using the equation y = mx + b
        y = slope * x + intercept
        # Plot the line
        plt.plot(x, y, label=f'y = {slope}x + {intercept}')
    
    # Add labels and a legend
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Plot of Lines with Different Slopes and Intercepts')
    plt.xlim(0, 30)
    plt.grid(True)