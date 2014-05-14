function ifGood = checkCoeffRoots(AR, MA, allowedValue, checkZeros)
% This function checks for causality and invertibility of the ARMA model
% under consideration; it returns 1 if ARMA model passed the check and 0 -
% otherwise

% checkZeros flag is set to 1 if the check for zeros roots in the AR and MA 
% polynomials is required, otherwise this flag is set to 0

% allowedValue sets the highest allowed bound for the roots

    roots_AR = roots([1; -AR(:)]);
    roots_MA = roots([1; MA(:)]);

    ifGood = 0;
    if ~(all(abs(roots_AR) < allowedValue) && all(abs(roots_MA) < allowedValue))
        return;
    end
    if checkZeros
        if (length(AR) > 1 && all(abs(roots_AR) == 0))
            return;
        elseif (length(MA) > 1 && all(abs(roots_MA) == 0))
            return;
        end
    end
    
    ifGood = 1;
end