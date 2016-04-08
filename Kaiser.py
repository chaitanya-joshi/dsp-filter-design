# Design of BPF and BSF by Kaiser Windowing
#
from numpy import pi, log10, zeros, sin
from math import floor, ceil
import scipy
import scipy.fftpack
from scipy.signal import kaiser
import pylab
from matplotlib.pyplot import plot, show
import matplotlib.pylab as plt
#
# Filter specifications.
#
M_f = 4 # Filter number
qm = floor(0.1*M_f)
rm = M_f - 10*qm
BLm = 4 + 0.9*qm + 2*rm
BHm = BLm + 10
print BLm, BHm
Wp1, Wp2 = BLm, BHm
# Transition Band
del_w = 2.0
# Stop band
Ws1 = Wp1 - del_w
Ws2 = Wp2 + del_w
# Tolerances
del1 = 0.015
del2 = 0.015
# w for Kaiser window
W1 = (Ws1 + Wp1)/2.0
W2 = (Wp2 + Ws2)/2.0
f_sampling = 100.0 # This frequency is in kHz (1000 cyles/sec)
#
# Normalizing specifications.
w1 = W1*(2*pi)/f_sampling
w2 = W2*(2*pi)/f_sampling
del_w = del_w*(2*pi)/f_sampling
A = -20.0*log10(min(del1,del2))
if (A>50.0):
	beta = 0.1102*(A-8.7)
elif (A>=21.0 and A<=50.0):
	beta = 0.5842*pow((A-21),0.4) + 0.07886*(A-21)
else:
	beta = 0.0
M = ceil((A - 8.0)/(2.285*del_w))
print 'del_w: ', del_w
print 'A: ', A
print 'M: ', M
window = kaiser(int(M+1), beta)
W = [w1,w2,pi]
print 'W: ', W
# G = [0.0,1.0,0.0,0.0] # For Band Pass Filter
G = [1.0,0.0,1.0,0.0] # For Band Stop Filter
fft_len = 512
h = zeros(fft_len)
for n in xrange(int(M+1)):
	for k in xrange(3):
		if (n==M/2):
			h[n] = h[n] + (G[k] - G[k+1])*W[k]/pi
		else:
			h[n] = h[n] + (G[k] - G[k+1])*sin(W[k]*(n - M/2))/(pi*(n - M/2))
	if (beta!=0.0):
		h[n] = h[n]*window[n]
freq_res = abs(scipy.fft(h))
freqs = scipy.fftpack.fftfreq(fft_len, 1.0/f_sampling) # Frequency axis will have kHz (1000 cycles/sec) unit.
plt.plot(freqs,freq_res)
# plt.plot(freqs,20*scipy.log10(freq_res))
plt.plot((W1,W1), (0,1), 'k-', color='green')
plt.plot((W2,W2), (0,1), 'k-', color='green')
plt.plot((0,50), (1-del1,1-del1), 'k-', color='green')
plt.plot((0,50), (1+del1,1+del1), 'k-', color='green')
plt.plot((0,50), (del2,del2), 'k-', color='black')
plt.show()
# pylab.plot(freq,inp_spec,'x--')