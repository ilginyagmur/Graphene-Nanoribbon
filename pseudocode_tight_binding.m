%% Tight-binding Matrix
eps_s = -8.868; eps_p = 0;
H_ss = -6.769; H_sp = -5.580; H_pp_sig = -5.037; H_pp_pi = -3.033; 
H_matrix = []; atom2_list = [];

for i = 1:size(ac_atoms_chosen, 1)
    atom1 = ac_atoms_chosen(1, :); 
    type1 = atom1(3);
    atom1_listed = atom1; atom1_listed(1,3) = 0; 
    first_nn_loop = []; second_nn_loop = [];
    index_first = find(ismember(first_nn, atom1_listed, 'rows'));
    index_second = find(ismember(second_nn, atom1_listed, 'rows')); 
    for k = index_first + 1:size(first_nn, 1)
        if ~isequal(atom1_listed, first_nn(k, :)) && first_nn(k ,3) ~= 0 
            first_nn_loop = [first_nn_loop; first_nn(k, :)];
        else
            break
        end
    end
    for k = index_second + 1:size(second_nn, 1)
        if ~isequal(atom1_listed, second_nn(k, :))  && second_nn(k ,3) ~= 0 
            second_nn_loop = [second_nn_loop; second_nn(k, :)];
        else
            break
        end
    end
    H_row = [];
    for j = 1:size(uc_atoms_chosen, 1)
        atom2 = uc_atoms_chosen(j, :);
        type2 = atom2(3);
        if atom1 == atom2 
            h_s = eps_s; h_px = eps_p; h_py = eps_p; h_pz = eps_p;
            H_partial0 = diag([h_s, h_px, h_py, h_pz]);
            H_partial1 = 10.*eye; %random variables for tunning
            H_partial2 = 20.*eye;
            if ~any(ismember(first_nn_loop, atom2, 'rows')) ...
                    && ~any(ismember(second_nn_loop, atom2, 'rows')) 
                H_partial = H_partial0;
            elseif any(ismember(first_nn_loop, atom2, 'rows')) ...
                    && ~any(ismember(second_nn_loop, atom2, 'rows'))  
                H_partial = H_partial0 + H_partial1;
            elseif ~any(ismember(first_nn_loop, atom2, 'rows')) ...
                    && any(ismember(second_nn_loop, atom2, 'rows'))
                H_partial = H_partial0 + H_partial2;
            elseif any(ismember(first_nn_loop, atom2, 'rows')) ... 
                    && any(ismember(second_nn_loop, atom2, 'rows'))
                H_partial = H_partial0 + H_partial1+ H_partial2;
            end 
        end
        H_row = [H_row, H_partial];
        atom2_list = [atom2_list; atom2];
    end
    H_matrix = [H_matrix; H_row];
end

