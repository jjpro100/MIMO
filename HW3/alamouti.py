from time import time
from numpy import sqrt
import random
import matplotlib.pyplot as plt

def coherent_det2(x):
    y = 0
    if x>=2:
        y = 3
    elif x<2 and x>=0:
        y = 1
    elif x<0 and x>=-2:
        y = -1
    else:
        y = -3
    return y

pam_map = [-3, -1, 1, 3]

t=time()

N = 50000
EbNodB_range = [0,5,10,15,25]
itr = len(EbNodB_range)
ber = [None]*itr

for n in range (0, itr): 
 
    EbNodB = EbNodB_range[n]   
    EbNo=10.0**(EbNodB/10.0)
    noise_std = 1/sqrt(2*EbNo)
    noise_mean = 0
    
    no_errors = 0
    for m in range (0, N):
        tx_symbol1 = random.choice(pam_map)
        tx_symbol2 = random.choice(pam_map)
                
        noise1 = (random.gauss(noise_mean, noise_std)+
                1j*random.gauss(noise_mean, noise_std))
        noise2 = (random.gauss(noise_mean, noise_std)+
                1j*random.gauss(noise_mean, noise_std))
        
        
        ch_coeff1 = (random.gauss(0,1/sqrt(2))+
                    1j*random.gauss(0,1/sqrt(2)))
        ch_coeff2 = (random.gauss(0,1/sqrt(2))+
                    1j*random.gauss(0,1/sqrt(2)))
        
        rx_symbol1 =  ((1/sqrt(2))*tx_symbol1*ch_coeff1+ 
                       (1/sqrt(2))*tx_symbol2*ch_coeff2 + noise1)
        rx_symbol2 = (-(1/sqrt(2))*tx_symbol2*ch_coeff1+ 
                       (1/sqrt(2))*tx_symbol1*ch_coeff2 + noise2)
        
        estimate1 = (ch_coeff1.conjugate()*rx_symbol1+
                     ch_coeff2*rx_symbol2.conjugate())
        estimate2 = (ch_coeff2.conjugate()*rx_symbol1-
                     ch_coeff1*rx_symbol2.conjugate())
        
        det_symbol1 = coherent_det2(estimate1.real)
        det_symbol2 = coherent_det2(estimate2.real)
        
        no_errors += 1*(tx_symbol1 != det_symbol1)+1*(tx_symbol2 != det_symbol2)  
          
    ber[n] = 1.0*no_errors/(2*N)
    print("EbNodB:", EbNodB)
    print("Numbder of errors:", no_errors)
    print("Error probability:", ber[n])
        
plt.plot(EbNodB_range, ber, 'bo-')
plt.axis([0, 10, 0.001, 0.125])
#plt.xscale('linear')
plt.yscale('log')
plt.xlabel('EbNo(dB)')
plt.ylabel('BER')
#plt.grid(True)
plt.title('BPSK Modulation - Alamouti Scheme')
plt.show()
print(time() - t, "seconds")