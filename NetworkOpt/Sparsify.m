function A_sparsified = Sparsify(A,cut_edge)
% This code sparsifies the adjacency matrix 'A' by removing the eges with
% weigth smaller than 'cut_edge', then the new matriz is renormalized to
% keep the sum of all entries to be the same. The code returns
% 'A_sparsified', the sparsified adjacency matrix
    original_sum = sum(A(:));
    
    % Sparsify the matrix
    A_sparsified = A;
    A_sparsified(abs(A_sparsified) < cut_edge) = 0;
    
    % Calculate the sum of the remaining non-zero entries after sparsification
    remaining_sum = sum(A_sparsified(:));
    
    % Check if there are any non-zero entries left
    if remaining_sum ~= 0
        % Rescale the remaining non-zero entries to match the original sum
        scaling_factor = original_sum / remaining_sum;
        A_sparsified = A_sparsified * scaling_factor;
    end
end