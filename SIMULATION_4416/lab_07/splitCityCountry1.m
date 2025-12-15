function [capital,Between,country]=splitCityCountry1(s)
    country=extractAfter(s, ',');
    capital=extractBefore(s, ',');
    Between=extractBetween(s,',',',');
end