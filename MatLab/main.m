%clear
close all
clc

sample_per_sym = 10; 

%test_signal = gmsk(:, 1)';

gmskmodulator=comm.GMSKModulator;
gmskmodulator.BandwidthTimeProduct=.3;
gmskmodulator.SamplesPerSymbol=sample_per_sym;
gmskmodulator.BitInput=true;
gmskmodulator.PulseLength=3;
gmskmodulator.SymbolPrehistory=[-1 1];
t = 0:1e-3:1e-3; 
s = cos(2*pi*200*t) + i * sin(2*pi*200*t); 
data = [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1]'; 
%data = (rand(24, 1) > 0.8); 
%data = ones(24, 1); 
test_signal = gmskmodulator(data)';

test_signal = awgn(test_signal, 10);

%test_signal = [s, test_signal];

freq = GMSK_demodulator(test_signal); 

figure 
plot(test_signal)

sensitive = (pi/2) / 9600; 
gain = 1 /sensitive; 

freq = freq * gain; 
preamble_template = test_signal; 

plot(freq)
grid on;
signal = zeros(1, 2 * length(test_signal)); 
noise = test_signal * 0; 
test_signal = [noise, test_signal, noise]; 
test_signal = test_signal * 5000; 
index = 1; 
for i = 1:length(test_signal)
    signal(1, index) = real(test_signal(1, i)); 
    index = index + 1; 
    signal(1, index) = imag(test_signal(1, i));
    index = index + 1; 
end


fileID = fopen('~/dev_dsp/GMSK_demodulator/Includes/preamble.bin','w');
fwrite(fileID,signal,'int16');
fclose(fileID);

bitstream = zeros(1, fix(length(test_signal) / sample_per_sym) + 1);
index = 1;

zeros_signal_0 = [preamble_template, zeros(1,length(test_signal) - 1)]; 
zeros_signal_1 = [test_signal, zeros(1,length(preamble_template) - 1)]; 
fft_0 = fft(zeros_signal_0); 
fft_1 = fft(zeros_signal_1); 

mul_fft = fft_1.*conj(fft_0); 

xcorr_real = ifft(mul_fft); 
figure 
title("Correlation with FFT"); 
plot(abs(xcorr_real))

hold on; 
plot(abs(xcorr(preamble_template, test_signal)))

% for i = 1 : sample_per_sym : length(freq)
%     sum = 0; 
%     last = freq(1, i); 
%     for j = i + 1 : i + sample_per_sym - 1
%         if last > freq(1, j)
%             sum = sum - 1; 
%         else
%             sum = sum + 1; 
%         end
%         last = freq(1, j); 
%     end
%     if sum < 0 
%         bitstream(1, index) = 1;
%     else
%         bitstream(1, index) = 0; 
%     end
%     index = index + 1; 
% end