h1 = [0.1162 0.2325 0.93 0.2325 0.1662];
h2 = [0.29 0.5 0.58 0.5 0.29];
Nc = 16;

lambda1 = fft([h1 zeros(1, Nc-length(h1))]);
lambda2 = fft([h2 zeros(1, Nc-length(h2))]);

bpsk_map = [-1 1];

figure(1)
plot(h1);
figure(2)
plot(h2);

Sim_number = 100000;

errors1 = 0;
errors2 = 0;
snr = 10;
for k = 1:Sim_number
    % modulation
    d_freq = [];
    for i = 1:Nc
        d_freq = [d_freq bpsk_map(randi([1 2]))];
    end

    % Putting d in time domain
    d_time = length(d_freq)*ifft(d_freq);
    % Add CP
    x = [d_time(Nc-(length(h1)-1):Nc) d_time];

    % Convolution with channel
    c1 = conv(h1,x);
    c2 = conv(h2,x);

%     sig_var1 = mean(abs(c1).^2);
%     sig_var2 = mean(abs(c2).^2);   
%     noise_var1 = sig_var1*10^(-snr/10);
%     noise_var2 = sig_var2*10^(-snr/10);
%     y1 = c1 + sqrt(noise_var1/2)*randn(1,length(c1)) + j*sqrt(noise_var1/2)*randn(1,length(c1));
%     y2 = c2 + sqrt(noise_var2/2)*randn(1,length(c2)) + j*sqrt(noise_var2/2)*randn(1,length(c2));
    

    % Generate and add noise
    noise_var = mean(abs(x).^2)*10^(-snr/10);
    y1 = c1 + sqrt(noise_var/2)*randn(1,length(c1)) + j*sqrt(noise_var/2)*randn(1,length(c1));
    y2 = c2 + sqrt(noise_var/2)*randn(1,length(c2)) + j*sqrt(noise_var/2)*randn(1,length(c2));    

    % Remove cp
    y1 = y1(length(y1)-Nc+1:end);
    y2 = y2(length(y2)-Nc+1:end);
    
    % Recover d
    d_freq_r1 = conj(lambda1).*fft(y1);
    d_freq_r2 = conj(lambda2).*fft(y2);

    dem1 = [];
    dem2 = [];

    % Coherent detection
    for i = 1:length(d_freq_r1)
        dem1 = [dem1 coherent_det(real(d_freq_r1(i)))];
        dem2 = [dem2 coherent_det(real(d_freq_r2(i)))];
    end

    errors1 = errors1 + sum(abs(dem1 - d_freq)./2);
    errors2 = errors2 + sum(abs(dem2 - d_freq)./2);
end
    
ser1 = errors1/Sim_number/Nc;
ser2 = errors2/Sim_number/Nc;

function y = coherent_det(x)
    y = 0;
    if(x>=0)
        y = 1;
    else
        y = -1;
    end
end