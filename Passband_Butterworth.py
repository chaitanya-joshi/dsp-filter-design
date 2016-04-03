# Design of Butterworth Pass Band Filter
#
from numpy import sqrt, log, zeros, sin, cos, pi
from math import ceil
from sympy import symbols, expand
s, z = symbols('s z') # Here, x is z^(-1)
#
### Band Pass to Low Pass transformation
def WL(w,w0,b):
	return (w**2 - w0**2)/(b*w)
#
### Passband Filter Specifications
# All frequenciea are in kHz
#
# Pass band
Wp1 = 12.0
Wp2 = 22.0
# Stop band
Ws1 = Wp1 - 2.0
Ws2 = Wp2 + 2.0
# Tolerances
del1 = 0.15
del2 = 0.15
#
# Parameters for the transformation WL = (W^2 - W0^2)/(B.W)
W0 = sqrt(Wp1*Wp2)
B = Wp2 - Wp1
transform = (s**2 + W0**2)/(B*s)
WLs1 = WL(Ws1,W0,B)
WLs2 = WL(Ws2,W0,B)
# Corresponding LPF (Low Pass Filter) specs
Wp = 1.0
Ws = min(abs(WLs1),abs(WLs2))
D1 = 1.0/((1.0 -del1)**2) - 1.0
D2 = 1.0/(del2**2) - 1.0
N = ceil(log(D2/D1)/(2.0*log(Ws/Wp)))
Wc = Ws*pow(D2,-(0.5/N)) # (Ws/Wc)^2N = D2
Wc = round(Wc, 2)
print 'N: ', N, '\t', 'Wc: ', Wc
DC = pow(Wc,N)
ki = zeros(N/2)
a1 = zeros(N/2)
a2 = zeros(N/2)
l = len(ki)
denominator = 1.0
for m in xrange(l):
	ki[m] = 2.0*Wc*sin((2.0*m + 1)*pi/(2.0*N))
	ki[m] = round(ki[m], 2)
	print 'ki[%d]: ' %m, ki[m]
	# DC = DC/(1.0 + ki[m] + Wc**2)
	# a1[m] = 2.0*(Wc**2 - 1.0)/(1.0 + ki[m] + Wc**2)
	# a2[m] = (1.0 - ki[m] + Wc**2)/(1.0 + ki[m] + Wc**2)
	denominator = expand(denominator*(s**2 + ki[m]*s + Wc**2))
# numerator = DC*expand(pow((1.0 + x),int(N)))
numerator = Wc**N
HcS = numerator/denominator
print 'HcS: ', HcS
bilinear = (1.0 - pow(z,-1))/(1.0 + pow(z,-1))
HcS_Bandpass = HcS.subs(s, transform)
print 'HcS_Bandpass: ', HcS_Bandpass
Hz = HcS_Bandpass.subs(s, bilinear)
print 'Hz: ', Hz