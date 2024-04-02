function [unit_cell_matrix, new_a1, new_a2] = find_unit_cell(c1, c2, chosen_atom, ribbon_positions, A_or_B_atoms)
%inputs = trial vectors, chosen atom, atom positions in the ribbon, A (or
%B)atom positions depending on the chosen atom
% outputs = plot of the unit cell
new_a1 = chosen_atom; new_a2 = chosen_atom;
found_atom_a1 = false; found_atom_a2 = false;
eps = 10^(-2);

for i = 1:size(ribbon_positions)
   trial_a1 = chosen_atom + c1 * i;
   M = sortrows(abs(A_or_B_atoms - ones(size(A_or_B_atoms,1), 1)*trial_a1));
   M1 = (M <= eps);
   if any(ismember(M1, [1 1], 'rows'))
        new_a1 = trial_a1; 
        found_atom_a1 = true;
        break;
    end
end

for i = 1:size(ribbon_positions)
   trial_a2 = chosen_atom + c2 * i;
   M = sortrows(abs(A_or_B_atoms - ones(size(A_or_B_atoms,1), 1)*trial_a2));
   M2 = (M <= eps);
   if any(ismember(M2, [1 1], 'rows'))
        new_a2 = trial_a2; 
        found_atom_a2 = true;
        break;
    end
end

if found_atom_a1
    plot([chosen_atom(1), new_a1(1)], [chosen_atom(2), new_a1(2)], 'r-', 'LineWidth', 2);
    hold on
end

if found_atom_a2
    plot([chosen_atom(1), new_a2(1)], [chosen_atom(2), new_a2(2)], 'r-', 'LineWidth', 2);
    hold on 
end

% if found_atom_a1 && found_atom_a2
%     plot([chosen_atom(1), new_a1(1)], [chosen_atom(2), new_a1(2)], 'r-', 'LineWidth', 2);
%     hold on
%     plot([chosen_atom(1), new_a2(1)], [chosen_atom(2), new_a2(2)], 'r-', 'LineWidth', 2);
%     plot([new_a1(1), new_a2(1), new_a2(1), new_a1(1), new_a1(1)], [new_a1(2), new_a1(2), new_a2(2), new_a2(2), new_a1(2)], 'r-', 'LineWidth', 2);
%     hold off
% end
% hold off

unit_cell = [];
x_min = min(new_a1(1), new_a2(1)); x_max = max(new_a1(1), new_a2(1));
y_min = min(new_a1(2), new_a2(2)); y_max = max(new_a1(2), new_a2(2));

for i = 1:size(ribbon_positions, 1)
    atom = ribbon_positions(i, :);
    if (atom(1) >= x_min - eps) && ...
       (atom(1) <= x_max + eps) && ...
       (atom(2) >= y_min - eps) && ...
       (atom(2) <= y_max + eps)
        unit_cell = [unit_cell atom];
    end
end

unit_cell_matrix = transpose(reshape(unit_cell, [2, size(unit_cell, 2)/2]));
end