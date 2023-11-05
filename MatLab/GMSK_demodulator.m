function [s] = GMSK_demodulator(in)
    re = real(in); 
    im = imag(in); 

    diff = volk_32fc_x2_multiply_conjugate_32fc_u_sse3(in(1, 2:length(in)), in(1, 1:length(in) - 1)); 

    mtr = atan(imag(diff)./real(diff)); 

    plot(mtr); 
    grid on


%     re_dt = finite_diffs(re); 
%     im_dt = finite_diffs(im); 
% 
%     I_mul = re_dt .* im;
%     Q_mul = im_dt .* re; 
% 
%     sum_env = -I_mul + Q_mul; 
% 
%     sum_IQ = re.^2 + im.^2; 
% 
%     freq = sum_env ./ sum_IQ; 

%     figure
%     plot(freq); 
%     title('Freq. detected')
%     grid on
%     hold on
% 
%     %figure 
%     plot(real(in)); 
%     title('Real part of input signal')
%     grid on 

    s = mtr; 

end 