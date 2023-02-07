import numpy as np
import matplotlib.pyplot as plt
import math 

n = 10000
delta = 0.05
var = 1
the_range = 10

indexes = []
lengths = []
bins = []
cdf = []
X = math.sqrt(var/2)*np.random.randn(n, 2).view(np.complex128)   # Generating data from n RVs
X_abs = abs(X) # Getting magnitude out of samples
idx_sum = 0

for i in np.arange(0, the_range, delta):
    temp = X_abs[(X_abs>i) & (X_abs<=(i+delta))] # Splitting into bins
    idx_sum += len(temp) # Summing samples less or equal to t
    lengths.append((1/len(X_abs))*len(temp)) # Getting data for empirical distribution
    indexes.append(idx_sum) # Saving indexes were each bin starts
    cdf.append((1/len(X_abs))*idx_sum) # Adding data for cdf
    bins.extend(temp)

plt.figure(0)
plt.plot(np.arange(0, the_range, delta), lengths) 
plt.title('PDF')
plt.xlabel('Range')
plt.ylabel('Likelihood')
plt.figure(1)
plt.plot(np.arange(0, the_range, delta), cdf)
plt.title('CDF')
plt.xlabel('Range')
plt.ylabel('Probability')
plt.show()

print(len(bins))
print(len(X_abs))
