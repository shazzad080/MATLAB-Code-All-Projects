function [mean_value,max_value,min_value]=blockStats(v)
    if ~isrow(v)||numel(v)~=100
        error('Input must be a 1x100 row vector.');
    end
    reshaped_matrix=reshape(v,5,20);
    mean_value=mean(reshaped_matrix)
    max_value=max(reshaped_matrix)
    min_value=min(reshaped_matrix)
end