function charactMSFE = getARMAmsfeFinSample(h, AR, MA, sigma2)
% Computes the characteristic error (MSFE) associated to an ARMA forecasting with
% horizon h

    Phi     = LagOp([1; -AR(:)]);
    Theta   = LagOp([1; MA(:)]);

    Psi = cell2mat(toCellArray(mrdivide_my(Theta, Phi, 'Degree', h, 'RelTol', 1e-20, 'Window',2000, 'AbsTol', 1e-20)));

    res = 0;

    for i = 1:min(h, length(Psi))
        res = res + Psi(i)^2;
    end

    charactMSFE = sigma2 * res;

end
