function out = figurate_number(n, varargin)
    a= ["Linear", "Triangular", "Square", "Pentagonal", "Hexagonal", "Heptagonal", "Octagonal"];
    if nargin == 1
        s= 3;
    elseif nargin == 2
        x = varargin{1}; % the input string is taken into variable x
        s = find(a==x) + 1;
    else
        out= "Too many input arguments";
        return
    end
    out = (((s-2)*n*(n-1))/2) + n;
end