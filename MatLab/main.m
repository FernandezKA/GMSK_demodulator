%clear
close all
clc

sample_per_sym = 10; 

%test_signal = gmsk(:, 1)';

gmskmodulator = comm.GMSKModulator('BitInput',true,'PulseLength',1, ...
                          'SamplesPerSymbol',sample_per_sym);

data = [1, 0 , 1, 1, 1, 0]'; 

test_signal = gmskmodulator(data)';

freq = GMSK_demodulator(test_signal); 

sensitive = (pi/2) / 9600; 
gain = 1 /sensitive; 

freq = freq * gain; 

plot(freq)
grid on;

bitstream = zeros(1, fix(length(test_signal) / sample_per_sym) + 1);
index = 1;

for i = sample_per_sym/2 : sample_per_sym : length(freq)
    if freq(1, i) < 0 
        bitstream(1, index) = 1;
    else
        bitstream(1, index) = 0; 
    end
    index = index + 1; 
end