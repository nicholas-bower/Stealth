import ROOT
xMin = 1000
xMax = 3500
StBkgs = []
StBkgs.append( ROOT.TF1("fSt+","1/763105876.464/TMath::Power(x/13000.,(6.05882+.46944))",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt0","1/243937094.587/TMath::Power(x/13000.,6.05882)",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt-","1/77988523.373/TMath::Power(x/13000.,(6.05882-.46944))",xMin,xMax) )
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
StBkgs[1].SetTitle("BKG Shape")
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
label.append("1/x^{[p_{1}+\sigma]}")
label.append("1/x^{p_{1}}")
label.append("1/x^{[p_{1}-\sigma]}")

leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
i = 0
for StBkg in StBkgs:
    leg.AddEntry(StBkg,label[i],"LP")
    i += 1
leg.SetBorderSize(0);
leg.Draw("SAME")
c.Update()
c.Print("CMSDAS_BKG_shape.png")
