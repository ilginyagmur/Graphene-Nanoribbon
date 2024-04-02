clear

%% Graphene Lattice
a = 1; %lattice constant

% Primitive vectors
a1 = [a*3/2, a*sqrt(3)/2];
a2 = [a*3/2, -a*sqrt(3)/2];

x_cells = 25;
y_cells = 25;

% Basis vectors
b1 = [0, 0];  
b2 = [2*a, 0]; 

atom_positions = [];
A_atoms = []; B_atoms = []; %2 types of atoms in graphene

pos1_initial = (-x_cells-1) * a1 + (-y_cells - 1) * a2 + b1;
pos2_initial = (-x_cells-1) * a1 + (-y_cells - 1) * a2 + b2;
pos3_initial = -(-x_cells-1) * a1 - (-y_cells - 1) * a2 + b1;
pos4_initial = -(-x_cells-1) * a1 - (-y_cells - 1) * a2 + b2;
atom_positions = [pos1_initial; pos2_initial; pos3_initial; pos4_initial];
A_atoms = [pos1_initial; pos3_initial];
B_atoms = [pos2_initial; pos4_initial];

for i = -x_cells:x_cells
    for j = -y_cells:y_cells
        % Calculate the displacement of the current unit cell
        r1 = (i-1) * a1 + (j-1) * a2; 
        pos1 = b1 + r1; pos2 = b2 + r1; pos3 = b1 - r1; pos4 = b2 - r1;
        if any(ismember(atom_positions, [pos1] , 'rows')) || any(ismember(atom_positions, [pos2] , 'rows')) || ... 
            any(ismember(atom_positions, [pos3] , 'rows')) || any(ismember(atom_positions, [pos4] , 'rows'))
            atom_positions = [atom_positions];
            A_atoms = [A_atoms];
            B_atoms = [B_atoms];
        else
            atom_positions = [atom_positions; pos1; pos2; pos3; pos4];
            A_atoms = [A_atoms; pos1; pos3];
            B_atoms = [B_atoms; pos2; pos4]; 
        end
    end
end


A_atoms = sortrows(A_atoms); B_atoms = sortrows(B_atoms);

figure(1)
hold on
axis equal
scatter(A_atoms(:,1), A_atoms(:,2), 50, 'filled');
scatter(B_atoms(:,1), B_atoms(:,2), 50, 'filled');
title('Graphene Lattice','FontSize',14,'FontWeight','bold');
xlabel('x/a','FontSize',12)
ylabel('y/a','FontSize',12)

%% ARMCHAIR

%% Nanoribbon Generation
W = 6; %nanoribbon width
m_ac = 0; %slope of the cutting line
ac_ribbon_positions = ribbon_generator(W, m_ac, atom_positions);

figure(2)
hold on
axis equal
scatter(ac_ribbon_positions(:,1), ac_ribbon_positions(:,2), 50, 'filled');
xlim([0, 2*W]);
title('Armchair Edges','FontSize',14,'FontWeight','bold');
xlabel('x/a','FontSize',12)
ylabel('y/a','FontSize',12)

%% Unit Cell (Rectangular) 
chosen_atom = ac_ribbon_positions(344, :);
c1 = [1, 0]; c2 = [0, 3*sqrt(3)]; %Expected angles of the new primitive vectors

if any(ismember(A_atoms, chosen_atom, 'rows'))
    A_or_B_atoms = A_atoms;
end

if any(ismember(B_atoms, chosen_atom, 'rows'))
    A_or_B_atoms = B_atoms; 
end

[ac_unit_cell, ac_vector1, ac_vector2] = find_unit_cell(c1, c2, chosen_atom, ac_ribbon_positions, A_or_B_atoms);
xlim([chosen_atom(1) - 5, chosen_atom(1) + 5])
hold off

%% Neighbours in the Unit Cell
distance1 = 1; %distance between nearest neighbours
distance2 = sqrt(3); %distance between second nearest neighbours
[ac_atoms, ac_atoms_chosen, first_nn_ac, second_nn_ac] = find_neighbours(ac_unit_cell, ac_vector1, ac_vector2, distance1, distance2);

% the code does not work as intended in the zigzag case, needs further
% modifications 

% %% Zigzag 
% eps = 0.1;
% W = 4 ; %adding a tolerance due to rounding errors 
% m_zz= 1/sqrt(3); %slope of the cutting line
% zz_ribbon_positions = ribbon_generator(W, m_zz, atom_positions);
% 
% figure(3)
% hold on
% axis equal
% xlim([0, 2*W]);
% ylim([0, 2*W]);
% scatter(zz_ribbon_positions(:,1), zz_ribbon_positions(:,2), 50, 'filled');
% title('Zigzag Edges','FontSize',14,'FontWeight','bold');
% xlabel('x/a','FontSize',12)
% ylabel('y/a','FontSize',12)
% 
% %% Unit Cell (Rectangular) 
% chosen_atom = zz_ribbon_positions(136, :);
% c1 = [-3/2, 3/2*sqrt(3)]; c2 = [3/2, sqrt(3)/2]; %Expected angles of the new primitive vectors
% 
% if any(ismember(A_atoms, chosen_atom, 'rows'))
%     A_or_B_atoms = A_atoms;
% end
% 
% if any(ismember(B_atoms, chosen_atom, 'rows'))
%     A_or_B_atoms = B_atoms; 
% end
% 
% [zz_unit_cell, zz_vector1, zz_vector2] = find_unit_cell(c1, c2, chosen_atom, zz_ribbon_positions, A_or_B_atoms);
% xlim([chosen_atom(1) - 5, chosen_atom(1) + 5])
% hold off
% 
% %% Neighbours in the Unit Cell
% distance1 = 1; %distance between nearest neighbours
% distance2 = sqrt(3); %distance between second nearest neighbours
% [atoms_zz, zz_atoms_chosen, first_nn_zz, second_nn_zz] = find_neighbours(zz_unit_cell, zz_vector1, zz_vector2, distance1, distance2);
















