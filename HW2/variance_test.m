noise_var = 8;
%y1 = sqrt(noise_var/2)*randn(1,10000) + j*sqrt(noise_var/2)*randn(1,10000);
Nc = 50;
cmp = [];

n = 456;
sum = 0;
    for k=0:99
        hl = sqrt(noise_var/2)*randn() + j*sqrt(noise_var/2)*randn();
        sum = [sum (1/sqrt(Nc))*exp(j*2*pi*(k*n/Nc))*hl];
    end


v = var(sum);