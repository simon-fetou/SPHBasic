import os
import numpy as np

def write_vtk(filename:str, positions:np.ndarray, velocities:np.ndarray,
               accelerations:np.ndarray, densities:np.ndarray, pressures:np.ndarray):
#def write_vtk(filename:str, positions:np.ndarray, velocities:np.ndarray
#              , densities:np.ndarray, pressures:np.ndarray):
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
        



def write_vtk_with_ghost(filename, fluid_pos, ghost_pos):
    """
    Write fluid and ghost particle positions to a VTK file.

    Parameters:
    - filename: str, name of the VTK file to save.
    - fluid_pos: np.ndarray of shape (N_fluid, 2), fluid particle positions.
    - ghost_pos: np.ndarray of shape (N_ghost, 2), ghost particle positions.
    """

    # Ensure the output directory exists
    os.makedirs('ghostVTK', exist_ok=True)
    
    # Full path to the file
    filepath = os.path.join('ghostVTK', filename )
    
    with open(filepath, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Fluid and Ghost Particles\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Combine all positions
        all_positions = np.vstack((fluid_pos, ghost_pos))
        num_particles = all_positions.shape[0]

        # Write points
        f.write(f"POINTS {num_particles} float\n")
        for pos in all_positions:
            f.write(f"{pos[0]} {pos[1]} 0.0\n")  # 2D positions (z=0)

        # Write cell data (optional)
        f.write("\nCELLS 0 0\n")
        f.write("\nCELL_TYPES 0\n")

        # Write point data to differentiate fluid and ghost particles
        f.write("\nPOINT_DATA {}\n".format(num_particles))
        f.write("SCALARS Particle_Type int 1\n")
        f.write("LOOKUP_TABLE default\n")
        f.write(" ".join(["0"] * len(fluid_pos)))  # 0 for fluid particles
        f.write(" ")
        f.write(" ".join(["1"] * len(ghost_pos)))  # 1 for ghost particles
        f.write("\n")
