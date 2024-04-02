function ribbon_positions = ribbon_generator(W, m, atom_positions)
%inputs = width of the nanoribbon, slope of the line the nanoribbon is cut through, positions of the atom in the lattice
%outputs = corresponding nanoribbon 
ribbon_positions = [];
for i = 1:size(atom_positions, 1)
    x_point = atom_positions(i,1); y_point = atom_positions(i,2);
    y1_line = m * x_point;
    y2_line = m * x_point  + W;
    if (y_point >= min(y1_line, y2_line)-0.1) && (y_point <= max(y1_line, y2_line))
        ribbon_positions = [ribbon_positions; atom_positions(i, :)];
    end
end
ribbon_positions = sortrows(ribbon_positions);
end

function ac_ribbon_positions = arm_chair_generator(W, atom_positions)
%input = width of the nanoribbon, positions of the atom in the lattice
%outputs = corresponding nanoribbon with armchair edges
ac_ribbon_positions = [];
for i = 1:size(atom_positions, 1)
    if atom_positions(i, 2) >= 0 && atom_positions(i, 2) <= W
        ac_ribbon_positions = [ac_ribbon_positions; atom_positions(i, :)];
    end
end
ac_ribbon_positions = sortrows(ac_ribbon_positions);
end
