function [capital, country] = splitCityCountry(s)
    parts   = strsplit(s, ',');
    capital = strtrim(parts{1});
    country = strtrim(parts{2});
end