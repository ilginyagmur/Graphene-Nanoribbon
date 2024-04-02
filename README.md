# Graphene-Nanoribbon
The code generates nanoribbons from graphene lattice, finds the unit cell atoms and first/second neighbouring atoms

# Functions 
ribbon_generator = cuts the graphene lattice with given parameters (width, cutting angle) to generate graphene nanoribbons (GNRs).
find_unit_cell = when given an initial position in the nanoribbon atom list obtained from ribbon_generator, finds the unit cell.
find_neighbours = given the unit cell atom positions obtained from find_unit_cell, finds the unique atoms in the unit cell and their first and second neighbouring atoms

# Main Script 
Creates the graphene lattice 
Generates armchair and zigzag nanoribbons and finds their neighbouring atom lists 
Please note that zigzag nanoribbons needs further tunning as the functions do not sometimes work as intended. 

# Future Works
The aim is to create tight-binding matrices from the neighbouring lists and then to find the electronic band structure of the GNRs. I have added a pseudocode, named pseudocode_tight_binding, to have an idea how to generate the Hamiltonian matrix elements.




