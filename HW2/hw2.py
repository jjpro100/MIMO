import numpy as np
import matplotlib.pyplot as plt
import math 
import random

def coherent_det(x):
    y = 0
    if x>=0:
        y = 1
    else:
        y = -1
    return y

SNRs = [0, 2, 4, 6, 8, 10, 12, 14, 16]
h1 = [0.1162, 0.2325, 0.93, 0.2325, 0.1662]
h2 = [0.29, 0.5, 0.58, 0.5, 0.29]
Nc = 16
CP_size = 7
lambda1 = np.fft.fft(h1, Nc)
lambda2 = np.fft.fft(h2, Nc)

plt.figure(1)
plt.stem(h1)
plt.figure(2)
plt.stem(h2)

bpsk_map = [-1, 1]
Sim_number = 100000
SERs1 = []
SERs2 = []

for snr in SNRs:
    ser1 =0.0
    ser2 = 0.0
    errors1 = 0.0
    errors2 = 0.0
    for k in range(Sim_number):
        
        # Modulation
        d_freq = []
        for i in range(Nc):
            d_freq.append(random.choice(bpsk_map))

        # Putting d in time domain
        d_time = np.fft.ifft(d_freq)
        # Add CP
        x = np.concatenate((d_time[-CP_size:], d_time), axis=None)

        # Convolution with channel
        c1 = np.convolve(h1,x)
        c2 = np.convolve(h2,x)

        # Generate and add noise
        noise_var1 = np.mean(np.square(abs(c1)))*(10**(-snr/10))
        noise_var2 = np.mean(np.square(abs(c2)))*(10**(-snr/10))
        w1 = np.squeeze(math.sqrt(noise_var1/2)*np.random.randn(len(c1), 2).view(np.complex128))
        w2 = np.squeeze(math.sqrt(noise_var2/2)*np.random.randn(len(c2), 2).view(np.complex128)) 

        y1 = c1 + w1
        y2 = c2 + w2 

        # Remove cp
        y1 = y1[CP_size:(CP_size+Nc)]
        y2 = y2[CP_size:(CP_size+Nc)]

        # Recover d
        d_freq_r1 = np.conjugate(lambda1)*np.fft.fft(y1)
        d_freq_r2 = np.conjugate(lambda2)*np.fft.fft(y2)
    
        dem1 = []
        dem2 = []
    
        # Coherent detection
        for i in range(len(d_freq_r1)):
            dem1.append(coherent_det(np.real(d_freq_r1[i])))
            dem2.append(coherent_det(np.real(d_freq_r2[i])))
        
        errors1 += (np.asarray(dem1) != np.asarray(d_freq)).sum()
        errors2 += (np.asarray(dem2) != np.asarray(d_freq)).sum()
    
    ser1 = errors1/(Sim_number*Nc)
    ser2 = errors2/(Sim_number*Nc)
    SERs1.append(ser1)
    SERs2.append(ser2)
    print(ser1, ser2)

plt.figure(0)
plt.semilogy(SNRs,SERs1, label = 'h1 = [0.1162 0.2325 0.93 0.2325 0.1662]')
plt.semilogy(SNRs,SERs2, label = 'h2 = [0.29 0.5 0.58 0.5 0.29]')
plt.legend()
plt.title("BER vs SNRdb")
plt.show()
