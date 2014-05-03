function [psiJac, piJac] = getJacobianPsiPi(AR, MA, Psi, Pi, P_max)

    p = length(AR);
    q = length(MA);
    r = max(p, q);

    AR_lag = LagOp([1; -AR(:)]);
    MA_lag = LagOp([1; MA(:)]);

    
    dpsi_theta = zeros(q, length(Psi));
    for i = 1:q
        M = zeros(length(Psi), 1);
        M(i + 1) = 1;
        Psi_lag = LagOp([M]);
        temp = cell2mat(toCellArray(mrdivide_my(Psi_lag, AR_lag, 'Degree', P_max, 'RelTol', 1e-20,'Window',2000, 'AbsTol', 1e-20)));
        
        differ_size = length(Psi)-length(temp);
        if differ_size > 0
            for j = 1:differ_size
                temp = [temp 0];
            end
        else
            temp = temp(1:length(Psi));
        end
        dpsi_theta(i, :) = temp(1:length(Psi));
     end

    dpsi_phi = zeros(p, length(Psi));
    for i = 1:p
        S = getShiftMatrix(length(Psi), i);
        temp1 = S * Psi';
        Psi_lag = LagOp([temp1(1:end)]);
        
        temp = cell2mat(toCellArray(mrdivide_my(Psi_lag, AR_lag, 'Degree', P_max, 'RelTol', 1e-20,'Window',2000, 'AbsTol', 1e-20)));
        
        differ_size = length(Psi)-length(temp);
        if differ_size > 0
            for j = 1:differ_size
                temp = [temp 0];
            end
        else
            temp = temp(1:length(Psi));
        end
        dpsi_phi(i, :) = temp(1:length(Psi));
    end

    psiJac = [dpsi_phi' dpsi_theta'];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AR_lag = LagOp([1; -AR(:)]);
    MA_lag = LagOp([1; MA(:)], 'Tolerance', 1e-30);
  
    
    dpi_phi = zeros(p, length(Pi));
    for i = 1:p
        M = zeros(length(Pi), 1);
        M(i + 1) = -1;
        Pi_lag = LagOp([M]);
        temp = cell2mat(toCellArray(mrdivide_my(Pi_lag, MA_lag, 'Degree', P_max, 'Window',2000, 'RelTol', 1e-20, 'AbsTol', 1e-20)));
        
        differ_size = length(Pi)-length(temp);
        if differ_size > 0
            for j = 1:differ_size
                temp = [temp 0];
            end
        else
            temp = temp(1:length(Pi));
        end
        dpi_phi(i, :) = temp(1:length(Pi));
     end

    dpi_theta = zeros(q, length(Pi));
    for i = 1:q
        S = getShiftMatrix(length(Pi), i);
        temp1 = -S * Pi';
        Pi_lag = LagOp([temp1(1:end)], 'Tolerance', 1e-30);
        
        temp = cell2mat(toCellArray(mrdivide_my(Pi_lag, MA_lag, 'Degree', P_max, 'Window', 2000, 'RelTol', 1e-20, 'AbsTol', 1e-20)));
        
        differ_size = length(Pi)-length(temp);
        if differ_size > 0
            for j = 1:differ_size
                temp = [temp 0];
            end
        else
            temp = temp(1:length(Pi));
        end
        dpi_theta(i, :) = temp(1:length(Pi));
    end

    piJac = [dpi_phi' dpi_theta'];
    
    
    
    
end



