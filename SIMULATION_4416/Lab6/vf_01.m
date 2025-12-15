function out = vf_01(varargin)
    v=varargin
    n=nargin
    x= cell2mat(varargin)
    out= max(x)
end