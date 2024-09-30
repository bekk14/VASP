from pymatgen.io.vasp import Vasprun

# Load the vasprun.xml file
vasprun = Vasprun("vasprun.xml")

# Get the last structure and total energy
last_structure = vasprun.final_structure
total_energy = vasprun.final_energy  # Total energy (in eV)
total_elements = len(last_structure)  # Total number of elements

# Get the cell vectors
cell_vectors = last_structure.lattice.matrix  # 3x3 matrix for cell vectors

# Prepare cell vector components
ax, ay, az = cell_vectors[0]
bx, by, bz = cell_vectors[1]
cx, cy, cz = cell_vectors[2]

# Get the forces from the last structure
forces = vasprun.ionic_steps[-1]['forces']
# Initialize velocities
try:
    # Try to get the velocities from the last ionic step
    velocities = vasprun.ionic_steps[-1]['velocities']
except KeyError:
    # If 'velocities' key is missing, set velocities to zero for each element
    velocities = [[0.0, 0.0, 0.0]] * len(last_structure)

# Prepare the output
output_lines = []
output_lines.append(f"{total_elements}\n")
output_lines.append(f"   time=    0.000 (fs)  Energy= {total_energy:.5f} (eV) Cell_Vectors= {ax:.6f} {ay:.6f} {az:.6f} {bx:.6f} {by:.6f} {bz:.6f} {cx:.6f} {cy:.6f} {cz:.6f}\n")

charge_file="charge.dat"
with open(charge_file, 'r') as f:
        charges = [line.strip() for line in f.readlines()]

for line in output_lines:
    print(line, end='')
# Print the results
print("Element Symbol  Coordinate_x  Coordinate_y  Coordinate_z  Force_x  Force_y  Force_z")
for i, site in enumerate(last_structure):
    element_symbol = site.species_string  # Get the element symbol
    coords = site.coords                  # Get the coordinates
    force = forces[i]                     # Get the force on the ion
    velocity = velocities[i]

    # Format the output for each ion
    print(f"{element_symbol}  {coords[0]:.6f}  {coords[1]:.6f}  {coords[2]:.6f}  {force[0]:.6f}  {force[1]:.6f}  {force[2]:.6f} {velocity[0]:.6f}  {velocity[1]:.6f}  {velocity[2]:.6f}  {charges[i]}  0.000 0.000 0.000")
    l=f"{element_symbol}  {coords[0]:.6f}  {coords[1]:.6f}  {coords[2]:.6f}  {force[0]:.6f}  {force[1]:.6f}  {force[2]:.6f} {velocity[0]:.6f}  {velocity[1]:.6f}  {velocity[2]:.6f} {charges[i]}  0.000 0.000 0.000 \n"
    output_lines.append(l)
# Print the output lines
output_file="POSCAR.md"
with open(output_file, 'w') as f:
        f.writelines(output_lines)

