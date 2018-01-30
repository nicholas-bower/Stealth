import numpy as np
from scipy.integrate import quad
params = []
p0_fits = []
p1_fits = []
p1err_fits = []
with open('fit_choice_params.dat','r') as f:
  for ijt in f:
    params = ijt.strip('\n').split(' ')
    p0_fits.append(np.float32(params[0]))
    p1_fits.append(np.float32(params[1]))
    p1err_fits.append(np.float32(params[2]))
integrals = []
for i in range (0,4):
    I = quad(lambda st: 1./((st/13000.)**(p1_fits[i]*np.log(st))), 1500, np.inf)
  #  I=quad(lambda x: 1, 0 ,2)
    integrals.append(I[0])
frac_dev = []
for i in range (0,4):
    frac_dev.append(100*(1-(integrals[i]/integrals[1])))
for i in range (0,4):
    print str (i +2) + "jets: Integral = " + str(integrals[i]) + " || Fractional Deviation (from n=2) = " + str(frac_dev[i]) +  "% "
    print str(p1_fits[i])
