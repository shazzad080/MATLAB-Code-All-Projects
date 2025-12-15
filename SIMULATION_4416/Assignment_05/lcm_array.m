function result=lcm_array(arr)
    if any(arr<=0)||any(mod(arr,1)~=0)
        error('All elements in the array must be positive integers.');
    end
    result=arr(1);
    for i=2:length(arr)
        result=(result*arr(i))/gcd(result,arr(i));
    end
end
