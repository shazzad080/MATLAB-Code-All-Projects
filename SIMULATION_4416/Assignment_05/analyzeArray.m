function [neg_count,zero_indices,pos_nums,stats,sq_nums, rearranged_array] = analyzeArray(a)
    neg_count=sum(a<0,'all')
    zero_indices.linear=find(a==0)';
    [row_idx, col_idx]=find(a==0);
    zero_indices.row=row_idx';
    zero_indices.col=col_idx';
    zero_indices
    pos_nums=a(a>0)'
    stats.mean=mean(a,'all');
    stats.std_dev=std(a,0,'all');
    stats
    is_square_num=(a>=0)&(round(sqrt(a)).^2==a);
    sq_nums=a(is_square_num)'
    evens=a(mod(a,2)==0 & a~=0);
    zeros_list=a(a==0);
    odds=a(mod(a,2)~=0);
    rearranged_array=reshape([evens(:);zeros_list(:);odds(:)], size(a))
end
