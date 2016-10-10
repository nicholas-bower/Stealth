#!/usr/bin/python

import os
import sys
import array
import numpy as np
import ROOT
import argparse
from scipy.stats import chisquare

parser = argparse.ArgumentParser(description='Provide file to process')
parser.add_argument('-i','--input', help='Input file name',required=True)
parser.add_argument('-l','--stmin', default=1000., help='min sT value',type=float)
parser.add_argument('-r','--stmax', default=3500., help='max sT value',type=float)
args = parser.parse_args()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add(args.input+".root")

nEntries = ggIn.GetEntries()
print "nEntries="+str(nEntries)
#nBins = 5
nBins = 20
xMin = args.stmin
xMax = args.stmax
nJtMin = 2
nJtMax = 7
# jet index for normalization/ratio comparison
iJtScale = 0 
# jet index used as baseline for determining bkg scaling of higher nJets
iJtBkg = 0
# bin index used as baseline for determining bkg scaling of higher nJets (starts at 1 since 0:underflow)
iBinBkgLo = 1
iBinBkgHi = 2
chi2Tot = 0.

## MAIN ##
def main():
	print "For sT range",xMin,"->",xMax
	hST = []
	for j in range(nJtMin,nJtMax+1):
		hST.append( ROOT.TH1F("h"+str(j)+"jet",str(j)+"jet",nBins,xMin,xMax) )

	for h in hST:
		h.Sumw2()

	for jEvt in range(nEntries):

		# Initialize event
		iEvt = ggIn.LoadTree(jEvt)
		if iEvt < 0:
			break
		nb = ggIn.GetEntry(jEvt)
		if nb <= 0:
			continue
		if jEvt % 10000 == 0:
			print "Processing entry " + str(jEvt)

		#nJets = ggIn.b_nJets
		#evtSt = ggIn.b_evtSt
		nJets = ggIn.b_nJet
		evtSt = ggIn.b_evtST

		# Fill St for each jet multiplicity
		for iJet in range(0,nJtMax-nJtMin+1):
			if nJets >= nJtMax-iJet:
				hST[nJtMax-nJtMin-iJet].Fill(evtSt)
				break

	for iJet in range(0,nJtMax-nJtMin+1):
		print str(iJet+2)+"-jet evts: "+str(hST[iJet].Integral())

	# Write out histos
	hFile = ROOT.TFile("hFile.root","RECREATE")
	for h in hST:
		h.Write()

	## ==== DRAW PLOTS ==== ##

	c1 = ROOT.TCanvas("c1","c1",600,600)
	c1.SetBorderSize(0);
	c1.SetFrameBorderMode(0)
	ROOT.gStyle.SetTitleBorderSize(0)
	ROOT.gStyle.SetOptStat(0)

	pUp = ROOT.TPad("upperPad", "upperPad",.005, .270, .995, .995)
	pDn = ROOT.TPad("lowerPad", "lowerPad",.005, .005, .995, .270)

	pUp.Draw()
	pDn.Draw()

	pUp.SetMargin(12.e-02,3.e-02,5.e-03,2.e-02)
	pDn.SetMargin(12.e-02,3.e-02,29.e-02,4.e-02)

	## Draw upper pad

	pUp.cd()
	pUp.SetTicky()

	ROOT.gPad.SetLogy()

	# Get normalization factors
	norm = []
	print "==============="
	for h in hST:
		norm_ = 0.
		# only use entries from first bin 
		for iBin in range (1,2):
			norm_ += h.GetBinContent(iBin)
			print "sT="+str( h.GetXaxis().GetBinLowEdge(iBin) )+" : N="+str(h.GetBinContent(iBin))
		norm.append(norm_)
		print "==============="

	# Get scale factor for background estimates
	nRefEntries = []
	print "==============="
	# Get number of entries in reference bin/s
	for h in hST:
		nRefEntries_ = 0
		# only use entries from ~first few bins
		for iBin in range (iBinBkgLo,iBinBkgHi+1):
			nRefEntries_ += h.GetBinContent(iBin)
		nRefEntries.append(nRefEntries_)
	# Rescale histograms
	iJt = 0
	jtScale = []
	fScale = open("ScaleFactor.txt", "w")
	for h in hST:
		if nRefEntries[iJt] > 0:
			jtScale.append(nRefEntries[iJtBkg] / nRefEntries[iJt])
			print "scale factor: "+str(nRefEntries[iJtBkg])+" / "+str(nRefEntries[iJt])+" = "+str(jtScale[-1])
			#fScale.write( "%d %f\n" % (iJt+2,jtScale[-1]) )
			fScale.write( "%f\n" % jtScale[-1] )
			for iBin in range (1,nBins+1):
				print "sT="+str( h.GetXaxis().GetBinLowEdge(iBin) )+" : N="+str(h.GetBinContent(iBin))+" -> "+str(h.GetBinContent(iBin)*jtScale[-1])
		print "==============="
		iJt += 1
	fScale.close()

	# Normalize distributions 
	maxSTs = []
	iJt = 0
	for h in hST:
		h.Scale(1./norm[iJt])
		maxSTs.append(h.GetMaximum())
		iJt += 1
		#h.Write()

	hFile.Close()
	'''
	print "==============="
	for h in hST:
		for iBin in range (1,nBins):
			print "sT="+str( h.GetXaxis().GetBinCenter(iBin) )+" : N="+str(h.GetBinContent(iBin))
		print "==============="
	'''

	#hST[0].GetYaxis().SetRangeUser(0.,1.4*max(maxSTs))
	hST[0].GetYaxis().SetRangeUser(2.e-04,1.2)
	hST[0].GetXaxis().SetTitle("S_{T} [GeV]")
	hST[0].GetXaxis().SetTitleOffset(1.1)
	#hST[0].SetTitle("1#gamma, Ht > 700 GeV")
	hST[0].SetTitle("")


	hST[0].SetLineColor(1)
	hST[0].Draw("E")
	#hST[0].Draw("")
	count = 0
	for h in hST:
		if count > 0:
			h.SetLineColor(count+1)
			h.Draw("SAME E")
			c1.Update()
			#hST[i].Draw("SAME")
		count += 1

	# Draw legend
	count = 0
	label = ''
	leg = ROOT.TLegend(0.75,0.65,0.86,0.92)
	for h in hST:
		if count < nJtMax - nJtMin:
			label = str(count+2)+" jets"
		else:
			label = "#geq"+str(count+2)+" jets"
		#label = str(count+2)+" jets"
		leg.AddEntry(h,label,"LP")
		count += 1
	#ROOT.gStyle.SetBorderSize(0);
	leg.SetBorderSize(0);
	leg.Draw()
	#Draw labels
	tex = ROOT.TLatex()
	tex.DrawLatexNDC(0.19,0.16,"N="+str(nEntries));
	tex.DrawLatexNDC(0.19,0.10,"1#gamma, Ht > 700 GeV");

	## Draw lower pad

	pDn.cd()
	pDn.SetTicky()
	pDn.SetGridy()

	fUnity = ROOT.TF1("fUnity","[0]",xMin,xMax)
	fUnity.SetParameter( 0,1. )
	fUnity.GetYaxis().SetRangeUser(0.2,2.8)

	fUnity.GetXaxis().SetTitle("S_{T} [GeV]")
	fUnity.GetXaxis().SetLabelSize(0.1)
	fUnity.GetXaxis().SetTitleSize(0.12)
	fUnity.GetXaxis().SetTickLength(0.1)
	fUnity.GetXaxis().SetTitleOffset(1.14)
	fUnity.GetYaxis().SetTitle("n_{j} / n_{j} = "+str(iJtScale+2))
	fUnity.GetYaxis().SetNdivisions(305)
	fUnity.GetYaxis().SetLabelSize(0.1)
	fUnity.GetYaxis().SetTitleSize(0.11)
	fUnity.GetYaxis().SetTickLength(0.04)
	fUnity.GetYaxis().SetTitleOffset(0.5)
	fUnity.SetLineColor(2)
	fUnity.SetLineWidth(1)
	fUnity.SetLineStyle(4)
	fUnity.SetTitle("")

	fUnity.Draw()

	gRatio = []
	# Fill ratio and errors
	iJt = 0
	#chi2Tot = 0
	for h in hST:
		if iJt == iJtScale:
			iJt += 1
			continue
		#gRatio.append( make_ratio_graph(str(iJt+2)+"jt_"+str(iJtScale+2)+"jt", h, hST[iJtScale], chi2Tot) )
		gRatio.append( make_ratio_graph("g"+str(iJt+2)+"jt_"+str(iJtScale+2)+"jt", h, hST[iJtScale]) )
		'''
		gRatio.append(ROOT.TGraphErrors())
		print "hST[" + str(iJt) + "]..." 
		ratioY = []
		for iBin in range (1,nBins+1):
			if not (hST[iJtScale].GetBinContent(iBin) > 0.):
				continue
			ratioX = h.GetXaxis().GetBinCenter(iBin)
			ratioY.append( h.GetBinContent(iBin)/hST[iJtScale].GetBinContent(iBin) )
			gRatio[-1].SetPoint( iBin, ratioX, ratioY[-1] )
			if not (h.GetBinContent(iBin) > 0.):
				continue
			errY = ratioY[-1]*np.sqrt( (h.GetBinError(iBin)/h.GetBinContent(iBin))**2 + (hST[iJtScale].GetBinError(iBin)/hST[iJtScale].GetBinContent(iBin))**2 )
			gRatio[-1].SetPointError( iBin, (xMax-xMin)/(2.*nBins), errY )
			#print " >>"+str(h.GetBinContent(iBin))+" "+str(hST[1].GetBinContent(iBin))
		chisq,pval = chisquare(ratioY,np.ones(len(ratioY)))
		#print str(iJt)+": "+str(chisquare(ratioY,np.ones(len(ratioY))).[0])
		print "chiSq["+str(iJt)+"]: "+str(chisq)
		'''
		gRatio[-1].SetMarkerStyle(20)
		gRatio[-1].SetMarkerColor(iJt+1)
		gRatio[-1].SetLineColor(iJt+1)
		gRatio[-1].SetLineWidth(1)
		gRatio[-1].SetMarkerSize(.9)
		gRatio[-1].Draw("P SAME")
		c1.Update()
		iJt += 1

	print "Total chi2 =",chi2Tot
	#c1.Print("sT.png")
	c1.Print("sT.eps")

#def make_ratio_graph(g_name, h_num, h_den, chi2Tot):
def make_ratio_graph(g_name, h_num, h_den):
	global chi2Tot
	print "Doing",g_name
	gae = ROOT.TGraphAsymmErrors()
	gae.SetName(g_name)
	h_rat = h_num.Clone()
	h_rat.Divide(h_den)
	for i in range(1, h_rat.GetNbinsX() + 1):
		n_rat = 0
		r_high = 0

		# tail = (1 - cl) / 2; for 95% CL, tail = (1 - 0.95) / 2 = 0.025
		tail = 0.16
		if h_num.GetBinError(i) == 0. or h_den.GetBinError(i) == 0.:
			continue
		n_num = pow(h_num.GetBinContent(i) / h_num.GetBinError(i), 2)
		n_den = pow(h_den.GetBinContent(i) / h_den.GetBinError(i), 2)
		q_low = ROOT.Math.fdistribution_quantile_c(1 - tail, n_num * 2,
				(n_den + 1) * 2)
		r_low = q_low * n_num / (n_den + 1)
		if n_den > 0:
			n_rat = n_num / n_den
			q_high = ROOT.Math.fdistribution_quantile_c(tail, (n_num + 1) * 2,
					n_den * 2)
			r_high = q_high * (n_num + 1) / n_den
			gae.SetPoint(i - 1, h_rat.GetBinCenter(i), h_rat.GetBinContent(i))
			gae.SetPointError(
					i - 1, h_rat.GetBinWidth(i) / 2, h_rat.GetBinWidth(i) / 2,
					n_rat - r_low, r_high - n_rat
					)
			gae.GetXaxis().SetRangeUser(h_rat.GetBinLowEdge(1),
					h_rat.GetBinLowEdge(h_rat.GetNbinsX()) +
					h_rat.GetBinWidth(h_rat.GetNbinsX()))
	fUnity = ROOT.TF1("fUnity","[0]",xMin,xMax)
	fUnity.SetParameter( 0,1. )
	gae.Fit("fUnity","M0","",xMin,xMax)
	print "chiSq / ndf =",fUnity.GetChisquare(),"/",fUnity.GetNDF()
	print "====================="
	chi2Tot += fUnity.GetChisquare()
	return gae

#_____ Call main() ______#
if __name__ == '__main__':
	main()
