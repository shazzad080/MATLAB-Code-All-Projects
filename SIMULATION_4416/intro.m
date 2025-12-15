clear all 
clc
x=input("Give the numbers: ");
if x<0
    disp("Invalid grade must be greater than zero")
elseif x>100
    disp("Invalid grade must be less than 100")
elseif (x>=90) && (x<=100)
    disp("Your grade is A+")
elseif (x>=80) && (x<=89)
    disp("Your grade is A")
elseif (x>=70) && (x<=79)
    disp("Your grade is A-")
elseif (x>=60) && (x<=69)
    disp("Your grade is B+")
else 
    disp("Fail!")
end