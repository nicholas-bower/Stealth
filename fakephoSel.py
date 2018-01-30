import os
import sys
import numpy.core.numeric as np
import argparse
import ROOT

# Register command line options
parser = argparse.ArgumentParser(description='Run STEALTH selection.')
parser.add_argument('-s','--sel', default='A', help='Selection scheme: A or B.',type=str)
parser.add_argument('-p','--real', default=1, help='Number of real photons',type=int)
parser.add_argument('-f','--fake', default=1, help='Number of Fake Photons',type=int)
parser.add_argument('-e','--era', default='', help='Run Era.',type=str)
parser.add_argument('-H','--ht',  default=60., help='HT cut.',type=float)
parser.add_argument('-r','--DeltaR',  default=0.4, help='DeltaR(pho,jet) cut.',type=float)
parser.add_argument('-n','--simNum',  default='', help='Simulation number', type =str)
parser.add_argument('-l','--iEvtMin', default = 0, help = 'starting Event number', type = int)
parser.add_argument('-m','--iEvtMax', default = 100000, help = 'Max Event number', type = int)
args = parser.parse_args()

# Initialize selection params
nPhoCut_  = 1
PhoEtLead = 20. 
PhoEtAll  = 20. 
if args.sel == 'A':
    nPhoCut_ = 2
    PhoEtLead = 35.
    PhoEtAll  = 25. 
    #PhoEtLead = 10. # MC
    #PhoEtAll  = 10. # MC
if args.sel == 'C':
    nPhoCut_ = 0
nFake    = args.fake
nReal    = args.real
HTcut_   = args.ht
runEra   = args.era
minDeltaRcut_ = args.DeltaR
iEventStart=args.iEvtMin
iEventEnd=args.iEvtMax
nJetsCut_ = 2
print " >> Running STEALTH 2016 Data Selection:",args.sel
print " >> HT cut:",HTcut_
print " >> Era: 2016%s"%runEra
print " >> DeltaR(pho,jt): %.2f"%minDeltaRcut_
## MAIN ##
def main():
#intitialize counters
##photon##              
                                                                                                                                                          
    njetsSel =0.
    cTotPhotons = 0.
    cRealSieie = 0.
    cFakeSieie = 0.
    cHoverE = 0.
    cNeuIso = 0.
    cREALCHIso = 0.
    cFAKECHIso = 0. 
    cFakePhoIso = 0.
    cRealPhoIso =0.
    cPhoID = 0.
    cPhoEta = 0.
    cNFake =0.
    cNReal =0.
    cSieieRealTest = 0.
    cSieieFakeTest = 0.
    cLeadEt=0.
    cPhoEle=0.
    cPhoEtAll=0.
    cPixelSeed = 0.

#incremental countsers
    iHoverE = 0.
    iNeuIso = 0.
    iFakeCHIso = 0.
    iRealCHIso = 0.
    iPhoIso = 0.
    iRealPhoIso =0.
    iPhoID = 0.
    iRealSieie = 0.
    iFakeSieie = 0.
    iRealLeadEt=0.
    iFakeLeadEt=0.
    iFakePhoEle=0.
    iRealPhoEle=0.
    iRealPhoEtAll=0.
    iFakePhoEtAll=0.
    iFakePixelSeed=0.
    iR9=0.
   # iPixelSeed = 0.
#jets                                                                                                                                                                                                                   
    cJets = 0.
    cJetEta = 0.
    cJetID = 0.
#Veto                
    cHLT = 0.
    cDRCut= 0.
    cFakeRej = 0.
    cRealRej =0.
    cNPhotons=0.
    cMu =0.
    cEle = 0.
    cPhoEta = 0.
    cNFRPhotons = 0.
    cBarrel = 0.
    cHT=0.
    rejReal = 0.
    rejFake=0.
    # Keep time
    sw = ROOT.TStopwatch()
    sw.Start()

    # For DeltaR cut
    minDeltaRij = 100.
    nJetsTot = 0.

    # Load input TTrees into TChain
   # eosDir = "/eos/cms/store/user/mandrews/"
   # ggInStr = "%s/DATA/ggSKIMS/DoubleEG_Run2016%s_ReminiAOD_HLTDiPho3018M90_SKIM*.root"%(eosDir,runEra)
    eosDir = "/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V08_00_26_04"
    ggInStr = "root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/skim-DoubleEG_FebReminiAOD.root"
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")
    #InputList = "/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH/inputFiles_2016%s.txt"%(runEra)
   # with open(InputList,'r') as f:
   #     for ijt in f:
   #         ggInStr = ijt.strip('\n')
   #         ggIn.Add(ggInStr)
    eosDir = "/eos/uscms/store/user/lpcsusystealth"
    ggIn.Add(ggInStr)
    nEvts = ggIn.GetEntries()
    print " >> Input file(s):",ggInStr
    print " >> nEvts:",nEvts

    # Initialize output file as empty clone
    eosOut = "/eos/user/n/nbower/Stealth/Ntuples/FakeReal/Jan20/2f"
    outFileStr = "%s/CutReminiAOD_Run2016%s_sel%s_evts_HT%d_evt%sto%s_%sf_%sr.root"%(eosOut,runEra,args.sel,HTcut_,iEventStart,iEventEnd,nFake,nReal) 
    outFile = ROOT.TFile(outFileStr, "RECREATE")
    outDir = outFile.mkdir("ggNtuplizer")
    outDir.cd()
    ggOut = ggIn.CloneTree(0) 
    print " >> Output file:",outFileStr

    # Initialize output branches 
    nJets_  = np.zeros(1, dtype=int)
    evtST_  = np.zeros(1, dtype=float)
    b_nJets = ggOut.Branch("b_nJets", nJets_, "b_nJets/I")
    b_evtST = ggOut.Branch("b_evtST", evtST_, "b_evtST/D")

    ##### EVENT SELECTION START #####

    # Event range to process
    iEvtStart = iEventStart
    iEvtEnd   = iEventEnd
 #   iEvtEnd   = 1000000
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

        evtST = 0.
        nTested+=1
        # Photon selection
        if ggIn.HLTPho>>14&1 == False: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
            cHLT+=1
            continue
        phoIdx  = [] # for DeltaR check: keep a list of photon indices passing photon selection
        nPhotons = 0.
        nFakepho = 0.
        chIso = 0.
        neuIso = 0.
        phoIso = 0.
        HoverE= 0.
        Sieie = 0.
        R9 = 0.
        passesFakeSieie = True
        passesFakechIso = True
        passesFakeEither = True
        LeadPho = True

############################Photon Selection###################################

        for i in range(ggIn.nPho):
            if abs(ggIn.phoEta[i]) > 1.442:
                cPhoEta +=1.
                continue
            cBarrel+=1
            HoverE = ggIn.phoHoverE[i]
            Sieie =  ggIn.phoSigmaIEtaIEtaFull5x5[i]
            R9 = ggIn.phoR9[i]
##########################Iso correction##############################                                                                                                                        
            if abs(ggIn.phoEta[i]) < 1.0:
                chIso = max(ggIn.phoPFChIso[i] - ggIn.rho*.0360, 0.)
                neuIso =  max(ggIn.phoPFNeuIso[i] - ggIn.rho*.0597, 0.)
                phoIso = max(ggIn.phoPFPhoIso[i] - ggIn.rho*.1210, 0.)
            else:
                chIso = max(ggIn.phoPFChIso[i] - ggIn.rho*.0377, 0.)
                neuIso =  max(ggIn.phoPFNeuIso[i] - ggIn.rho*.0807, 0.)
                phoIso = max(ggIn.phoPFPhoIso[i] - ggIn.rho*.1107, 0.)
#I'm trying a new implimentation of the photon specific curs because I'm not confident in my understanding of 
#python's implimentation of booleans so lets give this a whirl, its fairly similar to tanmay's implimentation
########################Fake specific cuts################################
            if Sieie > .01022 and Sieie < .015:
                passesFakeSieie = True
            else: 
                passesFakeSieie = False
            if (chIso > .441 and chIso < 15.):
                passesFakechIso = True
            else:
                passesFakechIso = False
            if passesFakechIso or passesFakeSieie:
                passesFakeEither=True
            else :
                passesFakeEither=False


####################Fake Photon Selection ##############                                                                                                                                                                
            if      (   ggIn.phoEt[i] > PhoEtAll
                    and neuIso < 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2
                    and phoIso < 2.571+0.0047*float(ggIn.phoEt[i])                                                                                     
                    and HoverE < .0396
                    and passesFakeEither== True
                    and ggIn.phohasPixelSeed[i] == 0  
                    and R9 < 1.0
                    ):
                if LeadPho == True and ggIn.phoEt[i] < PhoEtLead:
                    LeadPho = False
                    cLeadEt+=(ggIn.nPho-i)
                    break
                LeadPho = False
                nFakepho+=1
                cNFake+=1
                evtST += ggIn.phoEt[i]
                phoIdx.append(i)
            else: 
                rejFake+=1

####################Real Photon Selection ##############                                                                   
            if (ggIn.phoEt[i] > PhoEtAll
                    and ((ggIn.phoIDbit[i] >> 1) & 1) == 1                                                                    
                    and ggIn.phoEleVeto[i] == True                                                                                                                                                                                                   
                    ):
                if LeadPho == True and ggIn.phoEt[i] <PhoEtLead:
                    LeadPho = False
                    cLeadEt+=(ggIn.nPho-i)
                    break
                LeadPho=False
                nPhotons+=1
                cNReal+=1
                evtST += ggIn.phoEt[i]
                phoIdx.append(i)
            else:
                rejReal += 1
#############UpdateCounters######################
            cBarrel+=1
            if ggIn.phoEleVeto[i] == False:
                cPhoEle+=1
            if abs(ggIn.phoEta[i]) < 1.442:
                if ggIn.phoIDbit[i]>>1&1 == 1 and Sieie== .01022:
                    cSieieRealTest +=1
                if not ggIn.phoIDbit[i]>>1&1 == 1:
                    cPhoID+=1.
                    if Sieie== .01022:
                        cSieieFakeTest+=1
                if HoverE >.0306:
                    cHoverE+=1.
                if Sieie >.01022:
                    cRealSieie+=1.
                else:
                    cFakeSieie+=1.
                if neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2:
                    cNeuIso +=1.

          ##########Reversed counters#########

                if phoIso > 2.571+0.0047*float(ggIn.phoEt[i]):
                    cRealPhoIso+=1.
                if phoIso < 2.571+0.0047*float(ggIn.phoEt[i]) or phoIso<15:
                    cFakePhoIso += 1.
                if chIso <.441:
                    cREALCHIso+=1.
                if chIso<.441 or chIso> 15.:
                    cFAKECHIso +=1
                if ggIn.phoEt[i]> PhoEtAll:
                    cPhoEtAll+=1


#######################incremental counters######################
                if HoverE >.0306:
                    iHoverE +=1.
                if neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iNeuIso +=1.
                if phoIso > 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iPhoIso +=1.
                                           ############## Fake##########
                if passesFakeEither== False or phoIso > 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE> .0306:
                    iFakeSieie +=1.
                if ggIn.phoEt[i] < PhoEtAll or passesFakeEither== False or phoIso < 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iFakePhoEtAll +=1.
                if ggIn.phohasPixelSeed[i] != 0 or ggIn.phoEt[i] <= PhoEtAll or passesFakeEither== False or phoIso >= 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iFakePixelSeed += 1.
                if ggIn.phohasPixelSeed[i] != 0 or ggIn.phoEt[i] >= PhoEtAll or passesFakeEither== False or phoIso >= 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306 or R9>=1.0:
                    iR9 += 1.


                                           ####################REAL######################
                if Sieie > .01022 or phoIso > 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iRealSieie +=1.
                if (chIso > .441) or Sieie > .01022 or phoIso > 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iRealCHIso +=1
                if ggIn.phoEt[i] < PhoEtAll or (chIso > .441) or Sieie > .01022 or phoIso > 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE >.0306:
                    iRealPhoEtAll +=1.
                if ggIn.phoEleVeto[i] == False or ggIn.phoEt[i] < PhoEtAll or (chIso > .441) or Sieie > .01022 or phoIso > 2.571+0.0047*float(ggIn.phoEt[i]) or neuIso > 2.725+0.0148*float(ggIn.phoEt[i])+0.000017*float(ggIn.phoEt[i])**2 or HoverE>.0306:
                    iRealPhoEle += 1.



        if nFakepho!=nFake:  
            cFakeRej+=1
        if nReal!=nPhotons:
            cRealRej+=1
        if nPhotons != nReal or nFakepho != nFake:
            cNFRPhotons+=1

        if nPhotons != nReal or nFakepho != nFake:
            continue
        # Jet selection
        nJets = 0
        nJetsDR = 0
        evtHT = 0
        for i in range(ggIn.nJet):
            cJets+=1
            if (ggIn.jetPt[i] > 30.0
                    and abs(ggIn.jetEta[i]) < 2.4          
                    ):
                if not (ggIn.jetPFLooseId[i] != False and ggIn.jetPUID[i] > 0.61 and ggIn.jetID[i] == 6):
                    cJetID +=1
            if not (ggIn.jetPt[i] > 30.0
                    and abs(ggIn.jetEta[i]) < 2.4
                    and ggIn.jetPFLooseId[i] != False # for some reason == True doesnt work
                    and ggIn.jetPUID[i] > 0.62 # formerly jetPUidFullDiscriminant
                    and ggIn.jetID[i] == 6 # loose:2, tight:6
                    ):
                continue
            nJets += 1.
            evtHT += ggIn.jetPt[i] # Add jet pT to HT (even though not sure if it's photon)
            # DeltaR check: ensure this jet is well-separated from any of the good photons
            # To avoid double-counting, only add jet pT to ST if we're sure its not a photon 
            minDeltaRij = 100.
            for j in phoIdx: # loop over "good" photon indices
                dR = np.hypot(ggIn.phoEta[j]-ggIn.jetEta[i],ggIn.phoPhi[j]-ggIn.jetPhi[i]) #DeltaR(pho[j],jet[i])
                if dR < minDeltaRij: 
                    minDeltaRij = dR
            if minDeltaRij < minDeltaRcut_:
                continue
            nJetsDR += 1 # nJets passing the DeltaR check
            evtST += ggIn.jetPt[i]

        nJetsTot += nJetsDR
        if evtHT < HTcut_:
            cHT+=1
#            continue
        # Electron veto
        nEle = 0
        for i in range(ggIn.nEle):
            if not (ggIn.elePt[i] > 15.0
                    and ggIn.eleIDbit[i]>>3&1 == True # >>0:veto, >>1:loose, >>2:medium, >>3:tight
                    and abs(ggIn.eleEta[i]) < 2.5
                    and abs(ggIn.eleDz[i]) < 0.1
                    and ggIn.elePFPUIso[i] < 0.1
                    ):
                continue
            nEle += 1.
        if nEle != 0:
            cEle+=1.

        # Muon veto
        nMu = 0
        for i in range(ggIn.nMu):
            if not (ggIn.muPt[i] > 15.0
                    and ggIn.muPFPUIso[i] < 0.12
                    and ggIn.muIDbit[i]>>2&1 == True # >>0:loose, >>1:med, >>2:tight, >>3:soft, >>4:highpt
                    ):
                continue
            nMu += 1
        if nMu != 0:
            cMu+=1
        if nPhotons != nReal or nFakepho != nFake:
            continue
        if nMu != 0:
            continue
        if nEle != 0:
            continue
        # Photon selection                                                                                                                       
        # MET selection
        if ggIn.pfMET > 15.:
            evtST += ggIn.pfMET
        # Write this evt to output tree
        evtST_[0] = evtST
        nJets_[0] = nJetsDR
        ggOut.Fill()
        nAcc += 1.
        if nJetsDR >=2:
            njetsSel +=1

    ##### EVENT SELECTION END #####
    outFile.Write()
    outFile.Close()
   # cBarrel = cTotPhotons - cPhoEta-cLeadEt#number of photon candidates within the barrel with Et passing the leading cut)
    sw.Stop()
    totalFakeCut = cBarrel-cLeadEt-cNFake
    totalRealCut = cBarrel-cLeadEt-cNReal
    efHoverE = cHoverE/cBarrel*100
    efSieieReal = cRealSieie/cBarrel*100
    efNeuIso = cNeuIso/cBarrel*100
    efRealPhoIso = cRealPhoIso/cBarrel*100
    efFakeSieie = cFakeSieie/cBarrel*100
    with open('/afs/cern.ch/user/n/nbower/public/CMU/Stealth_SUSY/ST_plotting/CMSSW_8_0_24/src/STEALTH_PLAY/FakeReal/Selection_Outputs/Jan20/2f/Menglei2016%s_evt%sto%s_%sf_%sr.txt'%(runEra,iEvtStart,iEvtEnd,nFake,nReal),'w') as f:
        f.write(" Total Events in Era = " + str(nEvts)+"\n")
        f.write( " >> nAccepted evts:"+str(nAcc)+"/"+str(iEvtEnd-iEvtStart)+"("+str(100.*nAcc/(iEvtEnd-iEvtStart))+"% )"+"\n")
        f.write( " >> nJetsTot:"+str(nJetsTot)+"\n")
        f.write( " >> Real time:"+str(sw.RealTime()/60.)+"minutes"+"\n")
        f.write( " >> CPU time: "+str(sw.CpuTime() /60.)+"minutes"+"\n")
        f.write( "========================CounterResults==========================="+"\n")
        f.write( "####Photon Cut Effieciency#####"+"\n")
        f.write("LeadEt Cut = " + str(cLeadEt/(cBarrel)*100)+"%"+"\n")
        f.write("failing ID " +str(cPhoID)+"\n")
        f.write( "Total Photon candidates = " +str(cTotPhotons)+"\n")
       # f.write( "Barrel Efficiency = " +str(efBar) + "%"+"\n")
        f.write( "-Efficiency within Barel-"+"\n")
        f.write( "***REAL CUTS***"+"\n")
        f.write( "Photon ID = " + str(cPhoID/cBarrel*100)+"%"+"\n")
        f.write( "HoverE = " + str(efHoverE)+"%"+"\n")
        f.write( "SieieREAL = " + str(efSieieReal)+"%"+"\n")
        f.write( "Real Charged Hadron Iso = "+  str(cREALCHIso/cBarrel*100) + "%"+"\n")
        f.write( "Fake Charged Hadron Iso = "+  str(cFAKECHIso/cBarrel*100) + "%"+"\n")
        f.write( "Neutral Hadron Iso = " +str (efNeuIso)+ "%"+"\n")
        f.write( "Pho Iso = "+ str(efRealPhoIso)+"%"+"\n")
        f.write("Pho Et General cut (25 GeV) = " + str (cPhoEtAll/(cBarrel)*100)+"%"+"\n")
        f.write("Pho Electron Veto " + str(cPhoEle)+ " || " + str(cPhoEle/cBarrel*100)+ "%"+"\n") 
        f.write( "***Fake Cuts***"+"\n")
        f.write("SieieFAKE = " +  str(cFakeSieie/cBarrel*100)+"%"+"\n")
        f.write( "Pho Iso Fake = "+ str(cFakePhoIso/cBarrel*100)+"%"+"\n")
        f.write( "####Total number of Fake/ Real(passing efficiency)####"+"\n")
        f.write("Fake Photons = "+ str(cNFake)+"|| Passing Efficiency = " + str(cNFake/(cBarrel)*100)+"%"+"\n")
        f.write( "Real Photons = "+ str(cNReal)+"|| Passing Efficiency = " + str(cNReal/(cBarrel)*100)+"%"+"\n")
        f.write( "========================Incremental Counter Results==========================="+"\n")
        f.write( "-Efficiency within Barel-"+"\n")
        f.write( "***REAL CUTS***"+"\n")
        f.write( "HoverE = " + str(iHoverE/rejReal*100)+"%"+"\n")
        f.write( "Neutral Hadron Iso = " +str (iNeuIso/rejReal*100)+ "%"+"\n")
        f.write( "Pho Iso = "+ str(iPhoIso/rejReal*100)+"%"+"\n")
        f.write("Sieie = " +  str(iRealSieie/rejReal*100)+"%"+"\n")
        f.write("Charge Hadron Iso = " + str(iRealCHIso/rejReal*100) + "%"+"\n")
        f.write("Pho ET all = " + str (iRealPhoEtAll/rejReal*100) + "%"+"\n")
        f.write("Pho Electron Veto = " + str(iRealPhoEle/rejReal*100)+ "%"+"\n")
        f.write( "***FAKE CUTS***"+"\n")
        f.write( "HoverE = " + str(iHoverE/rejFake*100)+"%"+"\n")
        f.write( "Neutral Hadron Iso = " +str (iNeuIso/rejFake*100)+ "%"+"\n")
        f.write( "Pho Iso = "+ str(iPhoIso/rejFake*100)+"%"+"\n")
        f.write("Sieie = " +  str(iFakeSieie/rejFake*100)+"%"+"\n")
        f.write("Charge Hadron Iso = " + str(iFakeCHIso/rejFake*100) + "%"+"\n")
        f.write("Pho ET all = " + str (iFakePhoEtAll/rejFake*100) + "%"+"\n")
        f.write("Pho Pixel Veto = " + str(iFakePixelSeed/rejFake*100)+ "%"+"\n")
        f.write("R9 = " + str(iR9/rejFake*100)+ "%"+"\n")


        f.write( "===========Event Rejection=========="+"\n")
        f.write ( "HLT Eficiency=" + str(cHLT/nTested*100)+"%"+"\n")
        f.write( "Rejection due to combined number of fake and Real = "+ str(cNFRPhotons/(nTested-cHLT)*100)+"%"+"\n")
        f.write("HT<60GeV cut = " +str(cHT/(nTested-cNFRPhotons)*100)+"%" +"\n")
        f.write( "Muon Rejection = " + str(cMu/(nTested-cHLT-cNFRPhotons-cHT)*100)+"%"+"\n")
        f.write( "Electron Rejection = " + str (cEle/(nTested-cHLT-cHT-cNFRPhotons)*100)+"%"+"\n")
        f.write('=================Jet Selection=============\n')
        f.write('events passing jet selection = ' + str(njetsSel/nTested*100)+"%/n")
        f.write('Total Jet candidates '+ str(cJets)+'\n')
        f.write("Passing Cuts = " + str(cJetID/cJets*100)+"% \n") 
#_____ Call main() ______#
if __name__ == '__main__':
    main()
