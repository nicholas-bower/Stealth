import numpy as np
import ROOT
from scipy.integrate import quad
params = []
p0_fits = []
p1_fits = []
p1err_fits = []
with open('st_scaling_params.dat','r') as f:
  for ijt in f:
    params = ijt.strip('\n').split(' ')
    p0_fits.append(np.float32(params[0]))
    p1_fits.append(np.float32(params[1]))
    p1err_fits.append(np.float32(params[2]))
par_integrals = []
par_norms = []
for i in range (0,3):
  I = quad(lambda st: 1./((st/13000.)**((p1_fits[0]+(1-i)*p1err_fits[0])*np.log(st))), 1500, np.inf)
  N = quad(lambda st: 1./((st/13000.)**((p1_fits[0]+(1-i)*p1err_fits[0])*np.log(st))), 1300, 1500)
  par_integrals.append(I[0]/N[0])
  par_norms.append(N[0])
par_frac_dev = []
fit_frac_dev = []
normfit_frac_dev = []
N=[]
I=[]
I.append(quad(lambda st: 1.83648*10**-3/((st/13000.)**4.41145),1500, np.inf))
N.append(quad(lambda st:1.83648*10**-3/((st/13000.)**(4.41145)),1300, 1500))
I.append(quad(lambda st: 8.79246*10**-5/((st/13000.)**((p1_fits[0])*np.log(st))), 1500, np.inf))
N.append(quad(lambda st: 8.79246*10**-5/((st/13000.)**((p1_fits[0])*np.log(st))), 1300, 1500))
I.append(quad(lambda st:1.01805*10**3/np.exp(31.6012*st/13000.),1500, np.inf))
N.append(quad(lambda st:1.01805*10**3/np.exp(31.6012*st/13000.),1300, 1500))

fit_integrals = []
for i in range (0 , 3):
  fit_integrals.append(I[i][0]/N[i][0])
for i in range (0,3):
  par_frac_dev.append(100*((par_integrals[i]-par_integrals[1])/par_integrals[1]))
  normfit_frac_dev.append(100*(1-(fit_integrals[i]/fit_integrals[1])))
  fit_frac_dev.append(100*((I[i][0]-I[1][0])/I[1][0]))

print"===============Uncertainty due to Fit Parameter =================="
for i in range (0,3):
  print "Fractional deviation =" + str( par_frac_dev[i])+"%"
  print "N_"+str(i)+"="+str(par_norms[i])

print"===============Uncertainty due to normalized Fit Function==================="
for i in range (0,3):
  print "Fractional deviation =" + str( normfit_frac_dev[i])+"%"
  print "N_"+str(i)+"="+str(N[i][0])

print"===============Uncertainty due to Fit Function==================="
for i in range (0,3):
  print "Fractional deviation =" + str( fit_frac_dev[i])+"%"
  print "N_"+str(i)+"="+str(N[i][0])
xMin = 1300
xMax = 4000
StBkgs = []
StBkgs.append( ROOT.TF1("fSt+","(1/%d)/TMath::Power(x/13000.,((%d+%d)*TMath::Log(x)))"%(par_norms[0],p1_fits[0],p1err_fits[0]),xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt0","(1/%d)/TMath::Power(x/13000.,(%d*TMath::Log(x)))"%(par_norms[1],p1_fits[0]),xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt-","(1/%d)/TMath::Power(x/13000.,((%d-%d)*TMath::Log(x)))"%(par_norms[2],p1_fits[0],p1err_fits[0]),xMin,xMax) )
c = ROOT.TCanvas("c","c",600,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetLogy()
j = 0
lcolors = [4,1,2]
StBkgs[1].SetLineWidth(2)
StBkgs[1].SetLineColor(lcolors[1])
StBkgs[1].SetTitle("Background Estimation")
StBkgs[1].Draw()
for St in StBkgs:
    #St.SetParameter(iPar,float(St.GetParameter(iPar))/scales[i])                                                                                                                                                                                                              
  St.SetLineWidth(2)
  St.SetLineColor(lcolors[j])
  c.cd()
  St.Draw("SAME")
  c.Update()
  j += 1

label = []
label.append("1/x^{[p_{1}+\sigma]lnS_{t}}")
label.append("1/x^{p_{1}lnS_{t}}")
label.append("1/x^{[p_{1}-\sigma]lnS_{t}}")

leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
i = 0
for StBkg in StBkgs:
    leg.AddEntry(StBkg,label[i],"LP")
    i += 1
leg.SetBorderSize(0);
leg.Draw("SAME")
c.Update()
c.Print("sT_BKG_DiPho_Data_2016_std.png")
