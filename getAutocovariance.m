function autocov_vector = getAutocovariance(order, coeff_polynom, freq)
    autocov_vector = [];
    for i = 0:order
        s = 0;
        for j = 1:length(coeff_polynom) - i * freq
            s = s + coeff_polynom(j) * coeff_polynom(j + i * freq);
        end
        autocov_vector = [autocov_vector s]; 
    end
end

