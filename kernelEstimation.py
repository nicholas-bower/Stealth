import os
import sys
import numpy.core.numeric as np
import argparse
import ROOT
parser = argparse.ArgumentParser(description='Run STEALTH selection.')
parser.add_argument('-p','--real', default=1, help='Number of real photons',type=int)
parser.add_argument('-f','--fake', default=1, help='Number of Fake Photons',type=int)
parser.add_argument('-e','--era', default='', help='Run Era.',type=str)
parser.add_argument('-r', '--rho', default = 1, help = 'Rho for PDF', type=int)
args = parser.parse_args()
nFake    = args.fake
nReal    = args.real
runEra   = args.era
rho = args.rho
def main():

    # Load input TTrees into TChain
   
    eosDir = "/eos/user/n/nbower/Stealth/Ntuples/FakeReal/branched/"
    ggInStr = "%s/CutReminiAOD_Run2016*"%(eosDir)
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")
    ggIn.Add(ggInStr)
    
#    hFile = ROOT.TFile("ggInStr","READ")
    # Initialize output branches
    JetBranches = []
    b_Names = []
    b_Names.append("b_twoJetST")
    b_Names.append("b_threeJetST")
    b_Names.append("b_fourJetST")
    b_Names.append("b_fiveJetST")
    b_Names.append("b_sixpJetST")
    
    rooDataSets = []
    rooVars = []
    rooKeys = []
    mirror = ROOT.RooKeysPdf.MirrorLeft
#    binning = 
    for i in range(0,2):
        JetBranches.append(ROOT.gDirectory.Get(b_Names[i]))
        rooVars.append(ROOT.RooRealVar(b_Names[i],"rooVar",700,2000))
        rooDataSets.append(ROOT.RooDataSet(b_Names[i], b_Names[i], ggIn, ROOT.RooArgSet(rooVars[i])))
        rooKeys.append(ROOT.RooKeysPdf("keys"+ b_Names[i], "keys"+ b_Names[i], rooVars[i], rooDataSets[i], mirror, rho))
    c1 = ROOT.TCanvas("c1","c1",600,600)
    c1.SetBorderSize(0);
    c1.SetFrameBorderMode(0)
    xframe = rooVars[0].frame(700,2000,15)
    rooDataSets[0].plotOn(xframe)
   # xframe.Draw()
  #  c1.Print("Ktest")

    rooKeys[0].plotOn(xframe)
    xframe.Draw()
    c1.Print("Kernel_2Jet_1f1r_Left_Ro1.png")
    bframe = rooVars[1].frame(700,2000,15)

    rooDataSets[1].plotOn(bframe)
    #xframe.Draw()
    #c1.Print("Ktest")

    rooKeys[1].plotOn(bframe)
    xframe.Draw()
    c1.Print("Kernel3Jet_1f1r_Left_ro1.png")
    ##### EVENT SELECTION START #####
#_____ Call main() ______#
if __name__ == '__main__':
    main()
