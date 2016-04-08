# Design of Butterworth Pass Band Filter
#
from numpy import sqrt, log, zeros, sin, arccosh, pi, fft, polynomial, roots, complex64
import scipy
import scipy.fftpack
import pylab
from math import floor, ceil
from sympy import symbols, expand, simplify, fraction, Poly, cancel
# from sympy.plotting import plot
from matplotlib.pyplot import plot, show, errorbar
import matplotlib.pylab as plt
import sys
#
#
### Low Pass to Band Stop transformation
def WL(w,w0,b):
	return (b*w)/(w0**2 - w**2)
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
# Filter specifications.
M = 4 # Filter number
qm = floor(0.1*M)
rm = M - 10*qm
BLm = 4 + 0.9*qm + 2*rm
BHm = BLm + 10
print 'Stop band: ' '[', BLm, ',', BHm, ']'
s, x = symbols('s x') # Here, x is z^(-1)
#
f_sampling = 100.0 # This frequency is in kHz (1000 cyles/sec)
### Passband Filter Specifications
# All frequencies are in 1000 radians/sec
#
# Stop band
# Ws1 = 2.0*pi*BLm
# Ws2 = 2.0*pi*BHm
# # Pass band
# Wp1 = Ws1 - 2.0*pi*2.0
# Wp2 = Ws2 + 2.0*pi*2.0
# # Tolerances
# del1 = 0.15
# del2 = 0.15
#
### DUMMY SPECIFICATIONS FOR TESTING PURPOSES
# All frequencies are in 1000 radians/sec
#
# Stop band
Ws1 = 2.0*pi*35
Ws2 = 2.0*pi*40
# Pass band
Wp1 = Ws1 - 2.0*pi*2.0
Wp2 = Ws2 + 2.0*pi*2.0
# Tolerances
del1 = 0.15
del2 = 0.15
#
# Parameters for the transformation WL = (W^2 - W0^2)/(B.W)
W0 = sqrt(Wp1*Wp2)
B = Wp2 - Wp1
transform = (B*s)/(s**2 + W0**2)
# print transform
WLs1 = WL(Ws1,W0,B)
# print 'WLs1: ', WLs1
WLs2 = WL(Ws2,W0,B)
# print 'WLs2: ', WLs2
# Corresponding LPF (Low Pass Filter) specs (all frequencies are in 1000 radians/sec)
Wp = 1.0
Ws = min(abs(WLs1),abs(WLs2))
print 'Ws: ', Ws
D1 = 1.0/((1.0 -del1)**2) - 1.0
D2 = 1.0/(del2**2) - 1.0
# print 'D1: ', D1, 'D2: ', D2
# Evaluating epsilon and N for the Chebyshev LPF
eps = sqrt(D1)
print 'epsilon: ', eps
N = ceil(arccosh(sqrt(D2/D1))/(arccosh(Ws/Wp)))
print 'N: ', N
# Defining a list [0,0,...1] of size N+1 required to generate the coefficients of Nth order Cheyshev Polynomial
c = zeros(N+1)
c[-1] = 1.0
#
TnS = polynomial.chebyshev.Chebyshev(c)
TnS_coeff = polynomial.chebyshev.cheb2poly(TnS.coef)
# print 'Coefficients: ', TnS_coeff
CnW = 0.0
for i in xrange(int(N+1)):
	CnW = CnW + TnS_coeff[i]*(-1j*s)**i
# HcS_sq = 1.0/(1.0 + (eps**2)*(CnW**2))
denom_all = Poly(1.0 + (eps**2)*(CnW**2), s)
# print denom_all
denom_all = denom_all.all_coeffs()
leading_coeff = denom_all[0]
denom_all = [i/leading_coeff for i in denom_all]
norm_fac = (leading_coeff/abs(leading_coeff))*pow(abs(leading_coeff),0.5)
poles_all = roots(denom_all)
poles_LHP = zeros(N,dtype=complex64)
j = 0
denominator_LHP = 1.0 + 0j
for i in xrange(int(2*N)):
	if ((poles_all[i]).real<0):
		poles_LHP[j] = poles_all[i]
		j = j+1
		denominator_LHP = expand(denominator_LHP*(s-poles_all[i]))
# print poles
denom_poly = Poly(denominator_LHP,s)
denom_poly = denom_poly.all_coeffs()
denom_poly = [abs(i) for i in denom_poly]
denominator = 0.0
for i in xrange(len(denom_poly)):
	denominator = expand(denominator + denom_poly[int(N-i)]*s**(i))
# print poles_LHP
# print denom_poly
# denominator = simplify(denominator)
HcS = 1.0/(norm_fac*denominator)
# HcS = HcS.subs(s,1j*s)
HcS_Bandstop = HcS.subs(s,transform)
HcS_Bandstop = simplify(HcS_Bandstop)
# HcW_Bandstop = HcS_Bandstop.subs(s,1j*2.0*pi*s)
# plot(abs(HcW_Bandstop),(s,0,50))
bilinear = 2.0*f_sampling*(1.0 - x)/(1.0 + x)
Hz = HcS_Bandstop.subs(s, bilinear)
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
# ## Code to exit script. Used for debugging.
# sys.exit()
# ##