function output_matrix=flipped_diag(matrix,diag_type)
    n = size(matrix,1);
    output_matrix=matrix;
    if strcmp(diag_type,'main')
        main_diag_mask=logical(eye(n));
        main_diag_elements=output_matrix(main_diag_mask);
        output_matrix(main_diag_mask)=flip(main_diag_elements);
    elseif strcmp(diag_type,'anti')
        anti_diag_mask=logical(fliplr(eye(n)));
        anti_diag_elements=output_matrix(anti_diag_mask);
        output_matrix(anti_diag_mask)=flip(anti_diag_elements);
    end
end