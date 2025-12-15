function [add1, add2, add3] = joinword(s)
    add3 = extractAfter(s, ',') + ".com";
    add1 = extractBefore(s, ',') + ".com";
    add2 = extractBetween(s, ',', ',') + ".com";
end