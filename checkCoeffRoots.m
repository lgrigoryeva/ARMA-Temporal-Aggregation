function [ifGood] = checkCoeffRoots(AR, MA, allowedValue, checkZeros)

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