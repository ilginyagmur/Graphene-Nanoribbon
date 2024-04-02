function [atoms, uc_atoms_chosen, first_nn, second_nn] = find_neighbours(unit_cell, a1, a2, distance1, distance2)
%% inputs = atoms in the unit cell (N by 2 matrix), the vectors determining the unit cell,
%distance between the nearest neighbours, distance between the second
%nearest neighbours
%% outputs = unit cell with properties (type of the atom, if there is Bloch condition). 
% atoms chosen in the  unit cell, 
% nearest neighbours list (type=0 for chosen atom, type = 1 for corners, type = 2 
% for edge atoms, type = 3... for bulk atoms), second nearest neighbours
% (same types with first n.n.)
%% Unit Cell Atom Types
atoms = []; corner_atoms = []; edge_atoms = []; bulk_atoms = [];
eps = 10^(-5);
edge_type_count = 2;
bulk_type_count = 3;
dt = 10^(-2);

for i = 1:size(unit_cell, 1)
    chosen_atom = unit_cell(i, :);
    x = chosen_atom(1); 
    y = chosen_atom(2);
   if (abs(x - a1(1)) < eps && abs(y - a1(2)) < eps) ...
    || (abs(x - a2(1)) < eps && abs(y - a2(2)) < eps) ...
    || (abs(x - a2(1)) < eps && abs(y - a1(2)) < eps) ...
    || (abs(x - a1(1)) < eps && abs(y - a2(2)) < eps)
        type = 1; 
        isBloch = 1;
        corner_atoms = [corner_atoms; chosen_atom, type, isBloch];
   elseif (abs(x - a2(1)) < eps && min(a1(2), a2(2)) < y && y < max(a1(2), a2(2))) ...
           || (abs(x - a1(1)) < eps && min(a1(2), a2(2)) < y && y < max(a1(2), a2(2))) ...
           || (abs(y - a2(2)) < eps && min(a1(1), a2(1)) < x && x < max(a1(1), a2(1))) ...
           || (abs(y - a1(2)) < eps && min(a1(1), a2(1)) < x && x < max(a1(1), a2(1)))
       type = edge_type_count;
       edge_type_count = edge_type_count + dt;
       isBloch = 0;
       if (abs(x - a2(1)) < eps && min(a1(2), a2(2)) < y && y < max(a1(2), a2(2))) ...
               || (abs(x - a1(1)) < eps && min(a1(2), a2(2)) < y && y < max(a1(2), a2(2)))
           isBloch = 1;
       end
       edge_atoms = [edge_atoms; chosen_atom, type, isBloch];
    else
        type = bulk_type_count; 
        bulk_type_count = bulk_type_count + dt;
        isBloch = 0;
        bulk_atoms = [bulk_atoms; chosen_atom, type, isBloch];
    end
    atoms = [atoms; chosen_atom, type, isBloch];
end

atoms = sortrows(atoms ,3);
num_corner = size(corner_atoms, 1)*1/4; 
num_edge = size(edge_atoms, 1)*1/2;
num_bulk = size(bulk_atoms, 1)*1;

%% Choice of Atoms
%Edges, Corners, and Bulk
a1_vector = [a1, 0, 0] - corner_atoms(1, :); a2_vector = [a2, 0, 0] - corner_atoms(1, :);
pairs_list = [];

for i = 1:size(edge_atoms, 1)
    atom1 = edge_atoms(i, :);
    x1 = atom1(1); y1 = atom1(2);
    for j = i+1:size(edge_atoms, 1)
        atom2 = edge_atoms(j, :);
        x2 = atom2(1); y2 = atom2(2);
        if (abs(x2 - x1) < a1_vector(1) + eps && abs(x2 - x1) > a1_vector(1) - eps ...
           && abs(y2 - y1) < a1_vector(2) + eps && abs(y2 - y1) > a1_vector(2) - eps) ...
           || (abs(x2 - x1) < a2_vector(1) + eps && abs(x2 - x1) > a2_vector(1) - eps ...
           && abs(y2 - y1) < a2_vector(2) + eps && abs(y2 - y1) > a2_vector(2) - eps)
           pairs_list = [pairs_list; atom1; atom2];
        end
    end
end      

edge_choose = pairs_list(mod(1:size(pairs_list, 1), 2) == 1, :);
corner_choose = corner_atoms(1, :);
bulk_choose = bulk_atoms;
uc_atoms_chosen = [corner_choose; edge_choose; bulk_choose];

%% First Nearest Neighbours
first_nn = [];
first_vector1_list = []; 
    
for i = 1:size(uc_atoms_chosen, 1)
    atom1 = uc_atoms_chosen(i, :); type0 = 0;
    atom1_listed = atom1; 
    atom1_listed(1, 3) = type0;
    first_nn =[first_nn; atom1_listed]; 
    x1 = atom1(1); y1 = atom1(2);
    isBloch = atom1(4);
    if isBloch %if there is periodicity
        atom11 = atom1 + a1_vector; %the periodicity is always on a1 direction
        x11 = atom11(1); y11 = atom11(2);
    else
        x11 = x1; y11 = y1;
    end
    for j = 1:size(atoms, 1)
        atom2 = atoms(j, :);
        x2 = atom2(1); y2 = atom2(2);  
        vector1 = [x2 - x1, y2 - y1]; vector11 = [x2 - x11, y2 - y11];
        first_vector1_list = [first_vector1_list; vector1];
        min_range = vector11 - eps; max_range = vector11 + eps;
        d1 = sqrt((x1 - x2)^2 + (y1 - y2)^2); 
        d11 = sqrt((x11 - x2)^2 + (y11 - y2)^2);
        if abs(d1 - distance1) <= eps || ...
          (abs(d11 - distance1) <= eps && ...
         ~any(all(bsxfun(@ge, first_vector1_list, min_range) & ...
         bsxfun(@le, first_vector1_list, max_range), 2))) % If atom1 and atom11 gives the same type, exclude it
         first_nn = [first_nn; atom2];
        end
    end
end

%% Second Nearest Neighbours
second_nn = [];
second_vector1_list = []; 
    
for i = 1:size(uc_atoms_chosen, 1)
    atom1 = uc_atoms_chosen(i, :); type0 = 0;
    atom1_listed = atom1; 
    atom1_listed(1, 3) = type0;
    second_nn =[second_nn; atom1_listed];
    x1 = atom1(1); y1 = atom1(2);
    isBloch = atom1(4);
    if isBloch %if there is periodicity
        atom11 = atom1 + a1_vector; %the periodicity is always on a1 direction
        x11 = atom11(1); y11 = atom11(2);
    else
        x11 = x1; y11 = y1;
    end
    for j = 1:size(atoms, 1)
        atom2 = atoms(j, :);
        x2 = atom2(1); y2 = atom2(2);  
        vector1 = [x2 - x1, y2 - y1]; vector11 = [x2 - x11, y2 - y11];
        second_vector1_list = [second_vector1_list; vector1];
        min_range = vector11 - eps; max_range = vector11 + eps;
        d1 = sqrt((x1 - x2)^2 + (y1 - y2)^2); 
        d11 = sqrt((x11 - x2)^2 + (y11 - y2)^2);
        if abs(d1 - distance2) <= eps || ...
          (abs(d11 - distance2) <= eps && ...
         ~any(all(bsxfun(@ge, second_vector1_list, min_range) & ...
         bsxfun(@le, second_vector1_list, max_range), 2))) % If atom1 and atom11 gives the same type, exclude it
          second_nn = [second_nn; atom2];
        end
    end
end
end
