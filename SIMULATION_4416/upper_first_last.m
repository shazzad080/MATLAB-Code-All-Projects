function out=upper_first_last(s)
    out=lower(s);
    out(1)=upper(out(1));
    out(end)=upper(out(end));
end