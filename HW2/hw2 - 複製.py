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

def coherent_det2(x):
    out = []
    for e in x:
        if ((e-1).real**2 + e.imag**2) >= ((e+1).real**2 + e.imag**2):
            out.append(-1)
        else:
            out.append(1)
    out = np.array(out)
    return out

SNRs = [0, 2, 4, 6, 8, 10, 12, 14, 16]
h1 = [0.1162, 0.2325, 0.93, 0.2325, 0.1662]
h2 = [0.29, 0.5, 0.58, 0.5, 0.29]
Nc = 16

lambda1 = np.fft.fft(h1, Nc)
lambda2 = np.fft.fft(h2, Nc)

bpsk_map = [-1, 1]
Sim_number = 1000
SERs1 = []
SERs2 = []

errors1 = 0
errors2 = 0
snr = 10
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
        #print(d_freq)
        # Putting d in time domain
        d_time = np.fft.ifft(d_freq)
        # Add CP
        x = np.concatenate((d_time[-len(h1)-2:], d_time), axis=None)
        back = d_time[-7:]
        x = np.hstack((back,d_time))
        # Convolution with channel
        c1 = np.convolve(x,h1)
        c2 = np.convolve(h2,x)

        noise_var1 = np.mean((abs(c1**2)))*(10**(-snr/10))
        noise_var2 = np.mean(np.square(abs(c2)))*(10**(-snr/10))
        w1 = np.squeeze(math.sqrt(noise_var1/2)*np.random.randn(len(c1), 2).view(np.complex128))
        w2 = np.squeeze(math.sqrt(noise_var2/2)*np.random.randn(len(c2), 2).view(np.complex128))     
        w1 = math.sqrt(noise_var1/2)*(np.random.randn(len(c1))+1j*np.random.randn(len(c1)))

        # Generate and add noise
        # noise_var = np.mean(np.square(abs(x)))*(10**(-snr/10))
        # w1 = np.squeeze(math.sqrt(noise_var/2)*np.random.randn(len(c1), 2).view(np.complex128))
        # w2 = np.squeeze(math.sqrt(noise_var/2)*np.random.randn(len(c2), 2).view(np.complex128)) 

        y1 = c1 + w1
        y2 = c2 + w2 

        # Remove cp
        y1 = y1[7:(7+16)]
        
        y2 = y2[len(y2)-Nc:]

        # Recover d
        d_freq_r1 = np.conjugate(lambda1)*np.fft.fft(y1)
        
        d_freq_r2 = lambda2*np.fft.fft(y2)
    
        dem1 = []
        dem2 = []
    
        # Coherent detection
        # for i in range(len(d_freq_r1)):
        #     dem1.append(coherent_det(np.real(d_freq_r1[i])))
        #     dem2.append(coherent_det(np.real(d_freq_r2[i])))
        dem1 = coherent_det2(d_freq_r1)
        d_freq = np.asarray(d_freq)

        
        dem2 = coherent_det2(d_freq_r2)

        # errors1 += sum(abs(np.asarray(dem1) - np.asarray(d_freq))/2)
        # errors2 += sum(abs(np.asarray(dem2) - np.asarray(d_freq))/2)
        
        #errors1 += (np.asarray(dem1) != np.asarray(d_freq)).sum()
        for n in range(16):
            if dem1[n]!=d_freq[n]:
                errors1+=1        

        
        errors2 += (np.asarray(dem2) != np.asarray(d_freq)).sum()
    
    ser1 = errors1/(Sim_number*Nc)
    print(ser1)
    ser2 = errors2/(Sim_number*Nc)
    SERs1.append(ser1)
    SERs2.append(ser2)
    print(ser1, ser2)
plt.figure(0)
plt.semilogy(SNRs,SERs1)
plt.semilogy(SNRs,SERs2)
plt.show()
