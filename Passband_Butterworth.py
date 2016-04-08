# Design of Butterworth Pass Band Filter
#
from numpy import sqrt, log, zeros, sin, pi, fft
import scipy
import scipy.fftpack
import pylab
from math import floor, ceil
from sympy import symbols, expand, simplify, fraction, Poly, cancel
# from sympy.plotting import plot
from matplotlib.pyplot import plot, show, errorbar
import matplotlib.pylab as plt
import sys
s, x = symbols('s x') # Here, x is z^(-1)
#
# Filter specifications.
M = 4 # Filter number
qm = floor(0.1*M)
rm = M - 10*qm
BLm = 4 + 0.9*qm + 2*rm
BHm = BLm + 10
print BLm, BHm
#
### Low Pass to Band Pass transformation
def WL(w,w0,b):
	return (w**2 - w0**2)/(b*w)
#
def filter(ax, ay, x):
	m = len(ax)
	seq_len = len(x)
	y = zeros(seq_len)
	for i in xrange(seq_len):
		y[i] = ax[0]*x[i]
		fac = 1.0
		for j in xrange(m-1):
			if j>=i:
				fac = 0.0
			y[i] = y[i] + fac*(ax[j+1]*x[i-j-1] - ay[j+1]*y[i-j-1])
		# y[i] = (1.0/ay[0])*y[i]
	return y
#
f_sampling = 100.0 # This frequency is in kHz (1000 cyles/sec)
### Passband Filter Specifications
# All frequencies are in 1000 radians/sec
#
# Pass band
Wp1 = 2.0*pi*BLm
Wp2 = 2.0*pi*BHm
# Stop band
Ws1 = Wp1 - 2.0*pi*2.0
Ws2 = Wp2 + 2.0*pi*2.0
# Tolerances
del1 = 0.15
del2 = 0.15
#
### DUMMY SPECIFICATIONS FOR TESTING PURPOSES
# All frequencies are in 1000 radians/sec
#
# Pass band
# Wp1 = 2.0*pi*10.0
# Wp2 = 2.0*pi*40.0
# # Stop band
# Ws1 = Wp1 - 2.0*pi*2.0
# Ws2 = Wp2 + 2.0*pi*2.0
# # Tolerances
# del1 = 0.15
# del2 = 0.15
#
# Parameters for the transformation WL = (W^2 - W0^2)/(B.W)
W0 = sqrt(Wp1*Wp2)
B = Wp2 - Wp1
transform = (s**2 + W0**2)/(B*s)
# print transform
WLs1 = WL(Ws1,W0,B)
# print 'WLs1: ', WLs1
WLs2 = WL(Ws2,W0,B)
# print 'WLs2: ', WLs2
# Corresponding LPF (Low Pass Filter) specs (all frequencies are in 1000 radians/sec)
Wp = 1.0
Ws = min(abs(WLs1),abs(WLs2))
# print 'Ws: ', Ws
D1 = 1.0/((1.0 -del1)**2) - 1.0
D2 = 1.0/(del2**2) - 1.0
N = ceil(log(D2/D1)/(2.0*log(Ws/Wp)))
print 'N: ', N
Wc = Ws*pow(D2,-(0.5/N)) # (Ws/Wc)^2N = D2
# Wc = round(Wc, 2)
n_LHP = ceil(N/2)
ki = zeros(n_LHP)
a1 = zeros(n_LHP)
a2 = zeros(n_LHP)
denominator = 1.0
for m in xrange(int(n_LHP)):
	ki[m] = 2.0*Wc*sin((2.0*m + 1)*pi/(2.0*N))
	# ki[m] = round(ki[m], 2)
	# DC = DC/(1.0 + ki[m] + Wc**2)
	# a1[m] = 2.0*(Wc**2 - 1.0)/(1.0 + ki[m] + Wc**2)
	# a2[m] = (1.0 - ki[m] + Wc**2)/(1.0 + ki[m] + Wc**2)
	if (2*m+1 == N):
		denominator = expand(denominator*(s + Wc))
	else:
		denominator = expand(denominator*(s**2 + ki[m]*s + Wc**2))
 	# print denominator
# numerator = DC*expand(pow((1.0 + x),int(N)))
numerator = Wc**N
HcS = numerator/denominator
# freq_res_LPF = HcS.subs(s,1j*s)
# plot(abs(freq_res_LPF), (s,0,100))
# print 'HcS: ', HcS
HcS_Bandpass = HcS.subs(s, transform)
# freq_res_BPF = HcS_Bandpass.subs(s,1j*2.0*pi*s)
# plot(abs(freq_res_BPF), (s,0,50))
# print 'HcS_Bandpass: ', simplify(HcS_Bandpass)
bilinear = 2.0*f_sampling*(1.0 - x)/(1.0 + x)
Hz = HcS_Bandpass.subs(s, bilinear)
Hz = simplify(Hz)
Hz = cancel(Hz)
# print 'Hz: ', Hz
# print 'Numerator: ', fraction(Hz)[0]
# print 'Denominator: ', fraction(Hz)[1]
ax = Poly(fraction(Hz)[0], x) # Separating out the numerator
ax = ax.all_coeffs() # Array of size [(order of ax) + 1] with ax[0] being coeff of x^n
ax = list(reversed(ax)) # Reversing it so that ax[0] is coeff of order x^0
ay = Poly(fraction(Hz)[1], x) # Separating out the denominator
ay = ay.all_coeffs()
ay = list(reversed(ay))
norm_fac = ay[0]
print 'ay[0]: ', ay[0]
ax = [i/norm_fac for i in ax]
ay = [i/norm_fac for i in ay]
# print 'Length of ax: ', len(ax)
# print 'Length of ay: ', len(ay)
#
# ax[i] and ay[i] are the coefficients of x[n-i] and y[n-i]  respectively in the difference equation.
# The difference equation is implemented by the function 'filter' on any specified x[n]
# print 'Coefficients of numerator: ', ay
seq_len = 512 # 32768
h = zeros(seq_len)
inp_freq = 40.0 # This is in kHz (1000 cycles/sec)
h[0] = 1.0
# for i in xrange(seq_len):
# 	h[i] = sin(2.0*pi*inp_freq*i/f_sampling)
filtered_signal = filter(ax, ay, h)
# inp_spec = fft.fft(h)
inp_spec = abs(scipy.fft(h))
fft_len = len(abs(inp_spec))
print 'fft_len: ', fft_len
# freq = zeros(fft_len)
freq = scipy.fftpack.fftfreq(seq_len, 1.0/f_sampling) # Frequency axis will have kHz (1000 cycles/sec) unit.
# for i in xrange(fft_len):
# 	# freq[i] = i*(f_sampling/(2.0*pi))/seq_len
# 	freq[i] = (i-seq_len/2)*f_sampling/(seq_len)
# plt.plot(abs(inp_spec), color='b')
plt.plot(freq, abs(inp_spec), color='b')
# pylab.plot(freqs,20*scipy.log10(FFT),'x')
# pylab.plot(freq,inp_spec,'x--')
# plt.show()
# out_spec = fft.fft(filtered_signal)
out_spec = abs(scipy.fft(filtered_signal))
# plt.plot(abs(out_spec), color='r')
plt.plot(freq, abs(out_spec), color='r')
# pylab.plot(freq,out_spec,'x--')
# plt.plot((inp_freq,inp_freq), (0,200), 'k-', color='green')
plt.plot((Wp1/(2.0*pi),Wp1/(2.0*pi)), (0,1), 'k-', color='green')
plt.plot((Wp2/(2.0*pi),Wp2/(2.0*pi)), (0,1), 'k-', color='green')
plt.plot((Ws1/(2.0*pi),Ws1/(2.0*pi)), (0,1), 'k-', color='black')
plt.plot((Ws2/(2.0*pi),Ws2/(2.0*pi)), (0,1), 'k-', color='black')
plt.plot((0,50), (1-del1,1-del1), 'k-', color='green')
plt.plot((0,50), (del2,del2), 'k-', color='black')
plt.show()
# pylab.show()
## Code to exit script. Used for debugging.
sys.exit()
##