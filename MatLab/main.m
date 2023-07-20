%clear
close all
clc

sample_per_sym = 10; 

test_signal = gmsk(:, 1)'; 

freq = GMSK_demodulator(test_signal); 

freq_detect = zeros(1, 1); 

bit = 0; 


for i  = sample_per_sym : sample_per_sym : length(freq)
    if freq(1, i) > 0 
        bit = [bit, 1]; 
    else 
        bit = [bit, 0]; 
    end
end 

