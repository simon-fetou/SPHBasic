import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import re
import sys
# Adding the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Now you can import parameters as prm
import parameters as prm

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

def read_vtk(filepath):
    """
    Reads particle data (fluid and ghost) from ghost VTK file.
    
    Parameters:
    - filename: str, path to the VTK file.

    Returns:
    - fluid_positions: np.ndarray, positions of fluid particles.
    - ghost_positions: np.ndarray, positions of ghost particles.
    """
    fluid_positions = []
    ghost_positions = []
    particle_types = []

    with open(filepath, 'r') as f:
        lines = f.readlines()
        
        # Find where POINTS data starts
        for i, line in enumerate(lines):
            if line.startswith("POINTS"):
                num_points = int(line.split()[1])
                positions_start = i + 1
                break
        
        # Extract positions
        positions = []
        for j in range(positions_start, positions_start + num_points):
            pos = list(map(float, lines[j].strip().split()))
            positions.append(pos[:2])  # Only x, y for 2D
        
        # Find where POINT_DATA starts
        for i, line in enumerate(lines):
            if line.startswith("SCALARS Particle_Type"):
                particle_types_start = i + 2  # Data starts 2 lines below
                break
        
        
        # Extract particle types (may be on a single line or multiple lines)
        particle_data = lines[particle_types_start:particle_types_start + (num_points // 10) + 1]  # Approximate length
        particle_types = []

        for line in particle_data:
            values = line.strip().split()  # Split space-separated values
            particle_types.extend(map(int, values))  # Convert and append as integers

        # Ensure the length matches the number of points
        if len(particle_types) != num_points:
            raise ValueError(f"Mismatch: {len(particle_types)} particle types for {num_points} points.")



        # Separate positions into fluid and ghost particles
        for pos, ptype in zip(positions, particle_types):
            if ptype == 0:
                fluid_positions.append(pos)
            elif ptype == 1:
                ghost_positions.append(pos)

    return np.array(fluid_positions), np.array(ghost_positions)

def animate_particles(folder_path, domainBounds, save_as=None):
    """
    Animates fluid and ghost particles from a series of VTK files.

    Parameters:
    - vtk_folder: str, folder containing VTK files.
    - save_as: str (optional), filename to save the animation. If None, displays the animation.
    """

    xmin, xmax, ymin, ymax = domainBounds

    # Get all VTK files and sort them
    vtk_files = sorted(
                        [f for f in os.listdir(folder_path) if f.endswith('.vtk')],
                        key=lambda x: int(re.search(r'\d+', x).group())
                      )
    if not vtk_files:
        print("No VTK files found in the specified directory.")
        return
    
    # Initialize figure
    fig, ax = plt.subplots(figsize=(8, 8))
    fluid_scatter = ax.scatter([], [], c='blue', label='Fluid Particles', s=10)
    ghost_scatter = ax.scatter([], [], c='red', label='Ghost Particles', s=10,marker="x")

    # Domain boundaries
    ax.axvline(x=xmin, color='k', linestyle='--', linewidth=1.5)
    ax.axvline(x=xmax, color='k', linestyle='--', linewidth=1.5)
    ax.axhline(y=ymin, color='k', linestyle='--', linewidth=1.5)
    ax.axhline(y=ymax, color='k', linestyle='--', linewidth=1.5)
    ax.set_xlim(xmin - 0.2, xmax + 0.2)  # Adjust based on your domain size
    ax.set_ylim(ymin - 0.2, ymax + 0.2)  # Adjust based on your domain size
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(True)
    ax.set_title("Particle Visualization")
    ax.legend()


    def init():

        fluid_scatter.set_offsets(np.empty((0, 2)))  # Empty array with shape (0, 2) for 2D points
        ghost_scatter.set_offsets(np.empty((0, 2)))  # Same for ghost particles

        return fluid_scatter, ghost_scatter

    vtk_data = [read_vtk(os.path.join(folder_path, f)) for f in vtk_files]

    def update(frame):

        fluid_pos, ghost_pos = vtk_data[frame]
        fluid_scatter.set_offsets(fluid_pos)
        ghost_scatter.set_offsets(ghost_pos)
        ax.set_title(f"Particle Visualization - Time {frame * nsave * dt:.2f}")

        return fluid_scatter, ghost_scatter

    anim = FuncAnimation(fig, update, frames=len(vtk_files), init_func=init, repeat = False )#blit=False, interval=1000)
    
    plt.show()
    '''
    if save_as:
        ext = os.path.splitext(save_as)[1]
        writer = 'imagemagick' if ext == '.gif' else 'ffmpeg'
        anim.save(save_as, writer=writer, fps=10)
        print(f"Animation saved as {save_as}")
    '''

# Main application
if __name__ == "__main__":
    
    # Path to folder containing .vtk files
    vtk_folder_path = r"C:\Users\simon\GitProjects\SPHBasic\ghostVTK"

    domain_bounds=(x0D, lxDomain, y0D, lyDomain)
    anim = animate_particles(vtk_folder_path,domain_bounds, save_as="ghost_animation.mp4")
    
    plt.show()

