function [error] = getARMAmsfeFinSample(h, AR, MA, sigma2)

    Phi     = LagOp([1; -AR(:)]);
    Theta   = LagOp([1; MA(:)]);

    Psi = cell2mat(toCellArray(mrdivide_my(Theta, Phi, 'Degree', h, 'RelTol', 1e-20,'Window',2000, 'AbsTol', 1e-20)));

    res = 0;

    for i = 1:min(h, length(Psi))
        res = res + Psi(i)^2;
    end

    error = sigma2 * res;
%     error = 0;
end
