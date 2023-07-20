function [s] = GMSK_demodulator(in)
    re = real(in); 
    im = imag(in); 

    re_dt = finite_diffs(re); 
    im_dt = finite_diffs(im); 

    I_mul = re_dt .* im;
    Q_mul = im_dt .* re; 

    sum_env = -I_mul + Q_mul; 

    sum_IQ = re.^2 + im.^2; 

    freq = sum_env ./ sum_IQ; 

    figure
    plot(freq); 
    title('Freq. detected')
    grid on
    hold on

    %figure 
    plot(real(in)); 
    title('Real part of input signal')
    grid on 

    s = freq; 

end 