import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import sys
# Adding the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Now you can import parameters as prm
import parameters as prm
import re

#assigning global variables for local use

x0D = prm.x0D
y0D = prm.y0D
lxDomain = prm.lxDomain
lyDomain = prm.lyDomain
dx = prm.dx
dy = prm.dy

x0F = prm.x0F
y0F = prm.y0F
lxFluid = prm.lxFluid
lyFluid = prm.lyFluid

h = prm.h
alpha = prm.alpha

cRef = prm.cRef

m0 = prm.m0
k = prm.k

rho0 = prm.rho0

dt = prm.dt
Tf = prm.Tf
Nt = prm.Nt
nsave= prm.nsave

gamma = prm.gamma
g = prm.g
Npart = prm.Npart



def parse_vtk(filepath):
    """
    Function takes vtk files, classes its data in a list/dictionnary like

    Arguments:
    filepath: vtk local file path

    Returns:
    data: dictionnary of classified data
    """

    with open(filepath, 'r') as f:
        lines = f.readlines()               #Lists lines in the file
    
    #initializing a dict/list of data where will be stored each specific data
    data = {'points': [], 'velocity': [], 'acceleration':[], 'density': [], 'pressure': []}
    i = 0                                   # this counter will be used to go through lines
    num_points = 0                          # this counter will be used to go through 

    while i < len(lines):                       #go line by line adding +1 to the counter i
        line = lines[i].strip()                 # eliminates white spaces at the beg and end of line
        
        if line.startswith("POINTS "):
            num_points = int(line.split()[1])    # takes the total number of points stored
                                                 #Remember we wrote "POINTS numb" at the begin vtk file
            
            for j in range(i + 1, i + 1 + num_points):   #from next line where "POINTS" was found to last points line
                point_line = lines[j].strip().split()   

                if len(point_line) == 3:                 # Ensure valid 3D points (x,y,z)
                    try:
                        data['points'].append([float(v) for v in point_line]) #store (x,y,z) to points tag
                    except ValueError:
                        print(f"Skipping invalid point data: {lines[j]}")    #skip non_data lines if any
            i += num_points     # update the lines counter

        elif line.startswith("VECTORS Velocity"): #After points should start "POINTS_DATA" but when
                                                  # encounter "VECTORS Velocity" then ...

            # Same logic of storing as previously with POINTS
            for j in range(i + 1, i + 1 + num_points):
                velocity_line = lines[j].strip().split()

                if len(velocity_line) == 3:       # Ensure valid 3D vectors (Vx,Vy,Vz)
                    try:
                        data['velocity'].append([float(v) for v in velocity_line])
                    except ValueError:
                        print(f"Skipping invalid velocity data: {lines[j]}")
            i += num_points
            
        elif line.startswith("VECTORS Acceleration"):  #After Velocities should come Acceleration

            # Same logic as previously
            for j in range(i + 1, i + 1 + num_points):

                acceleration_line = lines[j].strip().split()
                if len(acceleration_line) == 3:  # Ensure valid 3D vectors (Ax,Ay,Az)
                    try:
                        data['acceleration'].append([float(v) for v in acceleration_line])
                    except ValueError:
                        print(f"Skipping invalid acceleration data: {lines[j]}")
            i += num_points
                
        elif line.startswith("SCALARS Density"):        #After Accel should start Densities
            i += 1                                      # Skip header and LOOKUP_TABLE line

            #Same logic as for vectors but these are SACALARS 
            for j in range(i+1, i+1 + num_points):
                density_line = lines[j].strip()
                try:
                    data['density'].append(float(density_line))
                except ValueError:
                    print(f"Skipping invalid density data: {lines[j]}")
            i += num_points

        elif line.startswith("SCALARS Pressure"):       #After Densities come Pressures
            i += 1                                      # Skip header and LOOKUP_TABLE line
            
            #Same logic as previously
            for j in range(i+1, i+1 + num_points):
                #print('Test:', lines[7212].strip())
                pressure_line = lines[j].strip()
                try:
                    data['pressure'].append(float(pressure_line))
                except ValueError:
                    print(f"Skipping invalid pressure data: {lines[j]}")
            i += num_points

        i += 1

    return data

def animate_vtk(folder_path):

    #vtk_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.vtk')])
    vtk_files = sorted(
                        [f for f in os.listdir(folder_path) if f.endswith('.vtk')],
                        key=lambda x: int(re.search(r'\d+', x).group())
                      )

    # Create subplots: One for density, one for pressure, and one for velocity magnitude
    fig = plt.figure(figsize=(12,8))
    plt.set_cmap('Blues')
    ax = fig.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    
    # Extract individual axes from the array
    ax_density = ax[0, 0]
    ax_pressure = ax[0, 1]
    ax_velocityQuiver = ax[1, 0]
    ax_velocity = ax[1, 1]

    # Initializing Global color limits
    densLimits = [float('inf'), -float('inf')]
    pressLimits = [float('inf'), -float('inf')]
    velLimits = [0, 0]

    # First pass to determine the global min/max values for each vtk
    for vtk_file in vtk_files:
        filepath = os.path.join(folder_path, vtk_file)
        data = parse_vtk(filepath)
        densLimits[0] = min(densLimits[0], min(data['density']))
        densLimits[1] = max(densLimits[1], max(data['density']))
        pressLimits[0] = min(pressLimits[0], min(data['pressure']))
        pressLimits[1] = max(pressLimits[1], max(data['pressure']))
        velMagnitude = np.linalg.norm(data['velocity'], axis=1)
        velLimits[0] = min(velLimits[0], min(velMagnitude))
        velLimits[1] = max(velLimits[1], max(velMagnitude))

    # Create colorbars for each plot
    density_cb = None
    pressure_cb = None
    velocity_cb = None
    acceleration_cb = None

    def update(frame):

        # Ensure these variables are accessed globally within the function
        nonlocal density_cb, pressure_cb, velocity_cb, acceleration_cb 

        # Clear the axes before each frame
        ax_density.clear()
        ax_pressure.clear()
        ax_velocityQuiver.clear()
        ax_velocity.clear()

        time_template = 'time = %.2fs'
        time_text = ax[0,0].text(1.1*lxDomain,1.2*lyDomain,'') 
        time_text.set_text(time_template%(0.))
        
        # Get the data for the current frame
        filepath = os.path.join(folder_path, vtk_files[frame])
        data = parse_vtk(filepath)

        # Extract the data from the parsed vtk data
        pos = np.array(data['points'])
        dens = np.array(data['density'])
        press = np.array(data['pressure'])
        vel = np.array(data['velocity'])
        #accel = np.array(data['acceleration'])
    
        # Compute the velocity magnitude
        velMagnitude =np.linalg.norm(vel, axis=1)

        #compute the acceleration magnitude
        #accel_magnitude = np.linalg.norm(accel, axis=1)
        
        # Plot the density
        scat1 = ax[0,0].scatter(pos[:,0],pos[:,1],c=dens,vmin=densLimits[0],vmax=densLimits[1],s=10)
        ax[0,0].grid(True)
        ax[0,0].set_xlim(x0D, lxDomain)
        ax[0,0].set_ylim(y0D, lxDomain)
        ax[0,0].set_xlabel('X')
        ax[0,0].set_ylabel('Y')

        if density_cb is None:  # Create the colorbar only once
            density_cb = fig.colorbar(scat1, ax=ax_density)
            density_cb.set_label('Density')

        # Plot the pressure
        scat2 = ax[0,1].scatter(pos[:,0],pos[:,1],c=press, vmin=pressLimits[0], vmax=pressLimits[1])
        ax[0,1].grid(True)
        ax[0,1].set_xlim(x0D, lxDomain)
        ax[0,1].set_ylim(y0D, lyDomain)
        ax[0,1].set_xlabel('X')
        ax[0,1].set_ylabel('Y')

        if pressure_cb is None:  # Create the colorbar only once
            pressure_cb = fig.colorbar(scat2, ax=ax_pressure)
            pressure_cb.set_label('Pressure')

        # Plot Velocity vectors
        quiv = ax[1,0].quiver(pos[:,0], pos[:,1], vel[:,0], vel[:,1], color=['r','b','g'], scale=21)
        ax[1,0].grid(True)
        ax[1,0].set_xlim(x0D, lxDomain)
        ax[1,0].set_ylim(y0D, lyDomain)
        ax[1,0].set_xlabel('X')
        ax[1,0].set_ylabel('Y')

        # Plot Velocity Magnitude
        scat3 = ax[1,1].scatter(pos[:,0], pos[:,1],c=velMagnitude,vmin=velLimits[0],vmax=velLimits[1],s=10)
        ax[1,1].grid(True)
        ax[1,1].set_xlim(x0D, lxDomain)
        ax[1,1].set_ylim(y0D, lyDomain)
        ax[1,1].set_xlabel('X')
        ax[1,1].set_ylabel('Y')

        if velocity_cb is None:  # Create the colorbar only once
            velocity_cb = fig.colorbar(scat3, ax=ax_velocity)
            velocity_cb.set_label('Velocity Magnitude')

        time_text.set_text(time_template%(frame*nsave*dt))

        return fig, scat1, scat2, quiv, scat3
    
    ani = animation.FuncAnimation(fig, update, frames=len(vtk_files), repeat=False)

    # Show the plot
    plt.show()


# Main application
if __name__ == "__main__":
    
    # Path to folder containing .vtk files
    vtk_folder_path = r"C:\Users\simon\GitProjects\SPHBasic\data"
    
    anim = animate_vtk(vtk_folder_path)
    plt.show()
