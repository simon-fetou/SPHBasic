import os
import numpy as np

def write_vtk(filename:str, positions:np.ndarray, velocities:np.ndarray,
               accelerations:np.ndarray, densities:np.ndarray, pressures:np.ndarray):
    
    """
    Function Writes SPH particle data to a .vtk file.
    
    Args:
        filename (str): Name of the output .vtk file.
        points (list of tuples): Particle positions [(x, y, z), ...].
        densities (list of float): Particle densities.
        pressures (list of float): Particle pressures.
        velocities (list of tuples): Particle velocities [(vx, vy, vz), ...].

    Returns:
            A .vtk file
    """

    # Ensure the output directory exists
    os.makedirs('data', exist_ok=True)
    
    # Full path to the file
    filepath = os.path.join('data', filename )

    with open(filepath, 'w') as f:
        # Header
        f.write("# vtk DataFile Version 4.2\n")
        f.write("SPH Particle Data\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Positions
        f.write(f"POINTS {len(positions)} float\n")

        for position in positions:
            f.write(f"{position[0]} {position[1]} {0.0}\n")    #Assuming 3D even if it's 2D
        
        # Position data
        f.write(f"POINTS_DATA {len(positions)}\n")
        
        # Velocity
        f.write("VECTORS Velocity float\n")
        for vel in velocities:
            f.write(f"{vel[0]} {vel[1]} {0.0}\n")            #Assuming 3D even if it's 2D
        
        # Acceleration
        f.write("VECTORS Acceleration float\n")
        for accel in accelerations:
            f.write(f"{accel[0]} {accel[1]} {0.0}\n")         #Assuming 3D even if it's 2D
        
        # Density
        f.write("SCALARS Density float\n")
        f.write("LOOKUP_TABLE default\n")
        for density in densities:
            f.write(f"{density}\n")
        
        # Pressure
        f.write("SCALARS Pressure float\n")
        f.write("LOOKUP_TABLE default\n")
        for pressure in pressures:
            f.write(f"{pressure}\n")
        
'''
# Example usage
points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0.5, 0.5, 0)]
densities = [1.0, 0.9, 1.2, 1.1, 1.0]
pressures = [1000.0, 950.0, 1200.0, 1100.0, 1050.0]
velocities = [(0.1, 0.2, 0.0), (0.0, -0.1, 0.0), (0.3, 0.4, 0.0), (-0.2, 0.0, 0.0), (0.0, 0.1, 0.0)]

write_vtk("sph_particles.vtk", points, densities, pressures, velocities)
'''