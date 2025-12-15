function result=vectorOps(v1,v2,operation)
    if ~isrow(v1)||~isrow(v2)
        error('Both inputs must be row vectors.');
    end
    if length(v1)~=length(v2)
        error('Input vectors must be of the same length.');
    end
    switch lower(operation)
        case 'add'
            result=v1+v2;
        case 'subtract'
            result=v1-v2;
        case 'multiply'
            result=v1.*v2;
        case 'divide'
            original_warning_state=warning('off', 'MATLAB:divideByZero');
            if any(v2==0)
                warning('Division by zero encountered. The result will contain NaN.');
            end
            result=v1./v2;
            result(v2==0)=NaN;
            warning(original_warning_state);
        otherwise
            error('Invalid operation. Supported operations are: add, subtract, multiply, divide.');
    end
end