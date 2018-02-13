import os
import sys
import numpy.core.numeric as np
import argparse
import ROOT
parser = argparse.ArgumentParser(description='Run STEALTH selection.')
parser.add_argument('-p','--real', default=1, help='Number of real photons',type=int)
parser.add_argument('-f','--fake', default=1, help='Number of Fake Photons',type=int)
parser.add_argument('-e','--era', default='', help='Run Era.',type=str)
parser.add_argument('-l','--iEvtMin', default = 0, help = 'starting Event number', type = int)
parser.add_argument('-m','--iEvtMax', default = -1, help = 'Max Event number', type = int)
args = parser.parse_args()
nFake    = args.fake
nReal    = args.real
runEra   = args.era
iEventStart=args.iEvtMin
iEventEnd=args.iEvtMax
def main():

    # Load input TTrees into TChain
   
    eosDir = "/eos/user/n/nbower/Stealth/Ntuples/FakeReal/Feb4/"
    ggInStr = "%s/CutReminiAOD_Run2016*"%(eosDir)
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")
    ggIn.Add(ggInStr)
    nEvts = ggIn.GetEntries()
    print " >> Input file(s):",ggInStr
    print " >> nEvts:",nEvts
    iEvtStart = iEventStart
    iEvtEnd   = iEventEnd
    if iEvtEnd == -1:
        iEvtEnd = nEvts
    # Initialize output file as empty clone
    eosOut = "/eos/user/n/nbower/Stealth/Ntuples/FakeReal/branched/"
    outFileStr = "%s/CutReminiAOD_Run2016%s_evts_HT60_evt%sto%s_%sf_%sr.root"%(eosOut,runEra,iEventStart,iEventEnd,nFake,nReal) 
    outFile = ROOT.TFile(outFileStr, "RECREATE")
    outDir = outFile.mkdir("ggNtuplizer")
    outDir.cd()
    ggOut = ggIn.CloneTree(0)
    nJets = 0
    evtST = 0.
    # Initialize output branches 
    twoJetST_  = np.zeros(1, dtype=float)
    threeJetST_  = np.zeros(1, dtype=float)
    fourJetST_ = np.zeros(1, dtype=float)
    fiveJetST_ = np.zeros(1, dtype=float)
    sixpJetST_ = np.zeros(1, dtype=float)
    b_twoJetST = ggOut.Branch("b_twoJetST", twoJetST_, "b_twoJetST/D")
    b_threeJetST = ggOut.Branch("b_threeJetST", threeJetST_, "b_threeJetST/D")
    b_fourJetST = ggOut.Branch("b_fourJetST", fourJetST_, "b_fourJetST/D")
    b_fiveJetST = ggOut.Branch("b_fiveJetST", fiveJetST_, "b_fiveJetST/D")
    b_sixpJetST = ggOut.Branch("b_sixpJetST", sixpJetST_, "b_sixpJetST/D")

    ##### EVENT SELECTION START #####

    # Event range to process
 #   iEvtStart = iEventStart
#    iEvtEnd   = iEventEnd
    print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"
    nTested= 0.
    nAcc = 0
    for jEvt in range(iEvtStart,iEvtEnd):

        # Initialize event
        if jEvt > nEvts:
            break
        treeStatus = ggIn.LoadTree(jEvt)
        if treeStatus < 0:
            break
        evtStatus = ggIn.GetEntry(jEvt)
        if evtStatus <= 0:
            continue
        if jEvt % 100000 == 0:
            print " .. Processing entry",jEvt
        

        nJets = ggIn.b_nJets
        evtST = ggIn.b_evtST
#sort st into branches by jets
        if nJets < 2:
            continue
        if nJets == 2:
            twoJetST_[0]=evtST
            ggOut.Fill()
            continue
        if nJets == 3:
            threeJetST_[0]=evtST
            ggOut.Fill()
            continue
        if nJets == 4:
            fourJetST_[0]=evtST
            ggOut.Fill()
            continue
        if nJets == 5:
            fiveJetST_[0]=evtST
            ggOut.Fill()
            continue
        if nJets >= 6:
            sixpJetST_[0]=evtST
            ggOut.Fill()
            continue
      #  ggOut.Fill()
    print(
    outFile.Write()
    outFile.Close()
#_____ Call main() ______#
if __name__ == '__main__':
    main()
