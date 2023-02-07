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


SNRs = [0, 5, 10, 15, 20, 25]

bpsk_map = [-1, 1]
pam_map = [-3, -1, 1, 3]
#pam_map = list((1/np.linalg.norm(pam_map))*np.asarray(pam_map))
Sim_number = 100
SERs_bpsk_m, SERs_bpsk_Al, SERs_pam_m, SERs_pam_Al = [], [], [], []
N_of_syms = 1200

def coherent_det4(x):
    i = np.argmin(abs(pam_map-x))
    return pam_map[i]


for snr in SNRs:
    err_bpsk_m = 0.0
    err_bpsk_Al = 0.0
    err_pam_m = 0.0
    err_pam_Al = 0.0
    for k in range(Sim_number):
        # Modulation
        bpsk_syms = []
        for i in range(N_of_syms):
            bpsk_syms.append(random.choice(bpsk_map))

        pam_syms = []
        for i in range(N_of_syms):
            pam_syms.append(random.choice(pam_map))

        for j in range(0, N_of_syms, 2):
            # Generate H matrix for Alamouti
            h1 = np.squeeze(math.sqrt(1/2)*np.random.randn(1, 2).view(np.complex128))
            h2 = np.squeeze(math.sqrt(1/2)*np.random.randn(1, 2).view(np.complex128))
            H = np.asarray([[h1, h2], [np.conj(h2), -np.conj(h1)]])
            H_conj = np.asarray(np.matrix(H).H)

            # Generate h1(t) and h2(t+1) and h1(t+2) and h2(t+3) for multi-shot, 
            # since we can just send one symbol in 2 slots when using 
            h1_m1 = np.squeeze(math.sqrt(1/2)*np.random.randn(1, 2).view(np.complex128))
            h2_m1 = np.squeeze(math.sqrt(1/2)*np.random.randn(1, 2).view(np.complex128))
            h1_m2 = np.squeeze(math.sqrt(1/2)*np.random.randn(1, 2).view(np.complex128))
            h2_m2 = np.squeeze(math.sqrt(1/2)*np.random.randn(1, 2).view(np.complex128))
            H_m1 = [[h1_m1, 0], [0, h2_m1]]
            H_m2 = [[h1_m2, 0], [0, h2_m2]]

            # Generate and add noise
            noise_var_bpsk = np.mean(np.square(bpsk_syms))*(10**(-snr/10))
            noise_var_pam = np.mean(np.square(pam_syms))*(10**(-snr/10))

            # BPSK noise
            w1 = np.squeeze(math.sqrt(noise_var_bpsk/2)*np.random.randn(1, 2).view(np.complex128))
            w2 = np.squeeze(math.sqrt(noise_var_bpsk/2)*np.random.randn(1, 2).view(np.complex128)) 

            w1_m1 = np.squeeze(math.sqrt(noise_var_bpsk/2)*np.random.randn(1, 2).view(np.complex128))
            w2_m1 = np.squeeze(math.sqrt(noise_var_bpsk/2)*np.random.randn(1, 2).view(np.complex128))

            w1_m2 = np.squeeze(math.sqrt(noise_var_bpsk/2)*np.random.randn(1, 2).view(np.complex128))
            w2_m2 = np.squeeze(math.sqrt(noise_var_bpsk/2)*np.random.randn(1, 2).view(np.complex128))  

            # Adding noise to bpsk symbols in multi-shot (slots 1 and 2)
            y1_m1 = h1_m1*bpsk_syms[j] + w1_m1
            y2_m1 = h2_m1*bpsk_syms[j] + w2_m1
            # (slots 3 and 4)
            y1_m2 = h1_m2*bpsk_syms[j+1] + w1_m2
            y2_m2 = h2_m2*bpsk_syms[j+1] + w2_m2

            # Adding noise to bpsk symbols in Alamouti
            y_Al = np.matmul(H, np.transpose([bpsk_syms[j], bpsk_syms[j+1]]))
            y_Al += np.transpose([w1, w2])

            # Equalization
            n1 = abs(h1_m1)*abs(h1_m1) + abs(h2_m1)*abs(h2_m1)
            n2 = abs(h1_m2)*abs(h1_m2) + abs(h2_m2)*abs(h2_m2)
            r1_m = (1/n1)*(np.conj(h1_m1)*y1_m1 + np.conj(h2_m1)*y2_m1)
            r2_m = (1/n2)*(np.conj(h1_m2)*y1_m2 + np.conj(h2_m2)*y2_m2)
            r_Al = (1/np.linalg.norm(H))*np.matmul(H_conj, y_Al)

            # Coherent detection
            if coherent_det(np.real(r1_m)) != bpsk_syms[j]: err_bpsk_m += 1
            if coherent_det(np.real(r2_m)) != bpsk_syms[j+1]:   err_bpsk_m += 1
            if coherent_det(np.real(r_Al[0])) != bpsk_syms[j]:  err_bpsk_Al += 1
            if coherent_det(np.real(r_Al[1])) != bpsk_syms[j+1]:   err_bpsk_Al += 1

            ########################################################################################################            

            # PAM noise
            EbNo=10.0**(snr/10.0)
            noise_std = 1/math.sqrt(2*EbNo)
            # noise_std = math.sqrt(noise_var_pam/2)
            w1 = np.squeeze(noise_std*np.random.randn(1, 2).view(np.complex128))
            w2 = np.squeeze(noise_std*np.random.randn(1, 2).view(np.complex128)) 

            w1_m1 = np.squeeze(noise_std*np.random.randn(1, 2).view(np.complex128))
            w2_m1 = np.squeeze(noise_std*np.random.randn(1, 2).view(np.complex128))

            w1_m2 = np.squeeze(noise_std*np.random.randn(1, 2).view(np.complex128))
            w2_m2 = np.squeeze(noise_std*np.random.randn(1, 2).view(np.complex128))  

            # Adding noise to bpsk symbols in multi-shot (slots 1 and 2)
            y1_m1 = h1_m1*pam_syms[j] + w1_m1
            y2_m1 = h2_m1*pam_syms[j] + w2_m1
            # (slots 3 and 4)
            y1_m2 = h1_m2*pam_syms[j+1] + w1_m2
            y2_m2 = h2_m2*pam_syms[j+1] + w2_m2

            # Adding noise to bpsk symbols in Alamouti
            y_Al = np.matmul(H, np.transpose([pam_syms[j], pam_syms[j+1]]))
            y_Al += np.transpose([w1, w2])

            # Equalization
            n1 = abs(h1_m1)*abs(h1_m1) + abs(h2_m1)*abs(h2_m1)
            n2 = abs(h1_m2)*abs(h1_m2) + abs(h2_m2)*abs(h2_m2)
            r1_m = (1/n1)*(np.conj(h1_m1)*y1_m1 + np.conj(h2_m1)*y2_m1)
            r2_m = (1/n2)*(np.conj(h1_m2)*y1_m2 + np.conj(h2_m2)*y2_m2)
            #r_Al = (1/np.linalg.norm(H))*np.matmul(H_conj, y_Al)
            r_Al = (1/(abs(h1)**2 + abs(h2)**2))*np.matmul(H_conj, y_Al)

            # Coherent detection
            if coherent_det4(np.real(r1_m)) != pam_syms[j]: err_pam_m += 1
            if coherent_det4(np.real(r2_m)) != pam_syms[j+1]:   err_pam_m += 1
            if coherent_det4(np.real(r_Al[0])) != pam_syms[j]:  err_pam_Al += 1
            if coherent_det4(np.real(r_Al[1])) != pam_syms[j+1]:   err_pam_Al += 1 
    
    ser_bpsk_m = err_bpsk_m/(Sim_number*N_of_syms)
    ser_bpsk_Al = err_bpsk_Al/(Sim_number*N_of_syms)
    ser_pam_m = err_pam_m/(Sim_number*N_of_syms)
    ser_pam_Al = err_pam_Al/(Sim_number*N_of_syms)

    SERs_bpsk_m.append(ser_bpsk_m)
    SERs_bpsk_Al.append(ser_bpsk_Al)
    SERs_pam_m.append(ser_pam_m)
    SERs_pam_Al.append(ser_pam_Al)

    print(snr)

plt.figure(0)
plt.semilogy(SNRs,SERs_bpsk_m, label = 'BPSK_multi_shot')
plt.semilogy(SNRs,SERs_bpsk_Al, label = 'BPSK_Alamouti')
plt.semilogy(SNRs,SERs_pam_m, label = '4PAM_multi_shot')
plt.semilogy(SNRs,SERs_pam_Al, label = '4PAM_Alamouti')
plt.legend()
plt.title("SER vs SNRdb")
plt.show()
