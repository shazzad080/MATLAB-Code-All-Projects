function category=number_category(n)
    possible_divisors=1:floor(n/2);
    proper_divisors=possible_divisors(mod(n,possible_divisors)==0);
    sum_of_divisors=sum(proper_divisors);
    if sum_of_divisors==n
        category='Perfect';
    elseif sum_of_divisors>n
        category='Abundant';
    else
        category='Deficient';
    end
end