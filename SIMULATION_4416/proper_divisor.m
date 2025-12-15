function out=proper_divisor(n)
    x=mod(n,1:n/2)==0
    out=find(x)
end