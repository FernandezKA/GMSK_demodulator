function cVector = volk_32fc_x2_multiply_conjugate_32fc_u_sse3(aVector, bVector)
    size = length(aVector); 
    cVector = zeros(1, size);
    for i = 1:length(aVector)
        cVector(1, i) = aVector(1, i) * conj(bVector(1, i)); 
    end 
end