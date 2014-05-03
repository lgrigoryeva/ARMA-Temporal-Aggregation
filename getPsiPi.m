function [Psi, Pi] = getPsiPi(AR, MA, P)
    AR_lag = LagOp([1; -AR(:)]);
    MA_lag = LagOp([1; MA(:)]);    
    
    Psi = cell2mat(toCellArray(mrdivide_my(MA_lag, AR_lag, 'Degree', P, 'RelTol', 1e-100, 'AbsTol', 1e-100, 'Window', 4000)));
    Pi = cell2mat(toCellArray(mrdivide_my(AR_lag, MA_lag, 'Degree', P, 'RelTol', 1e-100, 'AbsTol', 1e-100, 'Window', 4000)));
    for i = length(Psi) + 1:P + 1
        Psi(i) = 0;
    end
    for i = length(Pi) + 1:P + 1
        Pi(i) = 0;
    end
end