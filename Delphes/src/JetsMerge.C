//usage: root -l JetsMerge.C\(\"inputfile.root\"\,\"outputfile.root\"\)
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"
#include "classes/SortableObject.h"
#include "modules/Delphes.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTask.h"
#endif

#include <TVector.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>

void JetsMerge(const char *inputFile, const char *outputFile){
     gSystem->Load("libDelphes.so");
     TChain chain("Delphes");
     chain.Add(inputFile);
     ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
     TFile *file_sig = new TFile(inputFile);
     TFile *output = new TFile(outputFile, "recreate");
     TTree *tree_sig = (TTree*)file_sig->Get("Delphes");
     TTree *tree_output = new TTree("tree_output","Delphes");
     TLeaf *KTjet_size = tree_sig->GetLeaf("KTjet_size");
     TLeaf *GenJet_size = tree_sig->GetLeaf("GenJet_size");
     TLeaf *VLCjetR05N4_size = tree_sig->GetLeaf("VLCjetR05N4_size");

     TLeaf *VLCjetR05N4_eta = tree_sig->GetLeaf("VLCjetR05N4.Eta");
     TLeaf *VLCjetR05N4_phi = tree_sig->GetLeaf("VLCjetR05N4.Phi");
     TLeaf *VLCjetR05N4_pt = tree_sig->GetLeaf("VLCjetR05N4.PT");
     //TLeaf *VLCjetR05N4_btag = tree_sig->GetLeaf("VLCjetR05N4.BTag");
     TLeaf *VLCjetR05N4_deltaeta = tree_sig->GetLeaf("VLCjetR05N4.DeltaEta");
     TLeaf *VLCjetR05N4_deltaphi = tree_sig->GetLeaf("VLCjetR05N4.DeltaPhi");

     TLeaf *KTjet_eta = tree_sig->GetLeaf("KTjet.Eta");
     TLeaf *KTjet_phi = tree_sig->GetLeaf("KTjet.Phi");
     TLeaf *KTjet_pt= tree_sig->GetLeaf("KTjet.PT");
     //TLeaf *KTjet_btag = tree_sig->GetLeaf("KTjet.BTag");
     TLeaf *KTjet_deltaeta = tree_sig->GetLeaf("KTjet.DeltaEta");
     TLeaf *KTjet_deltaphi = tree_sig->GetLeaf("KTjet.DeltaPhi");

     TLeaf *GenJet_eta = tree_sig->GetLeaf("GenJet.Eta");
     TLeaf *GenJet_phi = tree_sig->GetLeaf("GenJet.Phi");
     TLeaf *GenJet_pt = tree_sig->GetLeaf("GenJet.PT");
     //TLeaf *GenJet_btag = tree_sig->GetLeaf("GenJet.BTag");
     TLeaf *GenJet_deltaeta = tree_sig->GetLeaf("GenJet.DeltaEta");
     TLeaf *GenJet_deltaphi = tree_sig->GetLeaf("GenJet.DeltaPhi");

     Int_t nEntries = tree_sig->GetEntries();
     
     Float_t KTGeneta[8];
     Float_t KTGenphi[8];
     Float_t KTGenpt[8];
     //Float_t KTGenbtag[8];
     Float_t KTGendeltaR[8];
     Int_t KTGenSize;
     
     Float_t VLCR05N4Genphi[4];
     Float_t VLCR05N4Geneta[4];
     Float_t VLCR05N4Genpt[4];
     //Float_t VLCR05N4Genbtag[4];
     Float_t VLCR05N4GendeltaR[4];
     Int_t VLCR05N4GenSize;

     tree_output->Branch("KTGenJetEta", &KTGeneta, "KTGeneta/F");
     tree_output->Branch("KTGenJetPhi", &KTGenphi, "KTGenphi/F");
     tree_output->Branch("KTGenJetPt", &KTGenpt, "KTGenpt/F");
     //tree_output->Branch("KTGenJetBtag", &KTGenbtag, "KTGenbtag/F");
     tree_output->Branch("KTGenJetDeltaR", &KTGendeltaR, "KTGenDeltaR/F");
     tree_output->Branch("KTGenJetSize", &KTGenSize, "KTGenSize/I");
     
     tree_output->Branch("VLCR05N4GenJetEta", &VLCR05N4Geneta, "VLCR05N4Geneta/F");
     tree_output->Branch("VLCR05N4GenJetPhi", &VLCR05N4Genphi, "VLCR05N4Genphi/F");
     tree_output->Branch("VLCR05N4GenJetPt", &VLCR05N4Genpt, "VLCR05N4Genpt/F");
     //tree_output->Branch("VLCR05N4GenJetBtag", &VLCR05N4Genbtag, "VLCR05N4Genbtag/F");
     tree_output->Branch("VLCR05N4GenJetDeltaR", &VLCR05N4GendeltaR, "VLCR05N4GendeltaR/F");
     tree_output->Branch("VLCR05N4GenJetSize", &VLCR05N4GenSize, "VLCR05N4GenSize/I");

     Float_t KTeta;
     Float_t KTphi;
     Float_t KTpt;
     //Float_t KTbtag;
     Float_t KTdeltaeta;
     Float_t KTdeltaphi;
     Float_t KTdeltaR;
     Float_t VLCeta;
     Float_t VLCphi;
     Float_t VLCpt;
     //Float_t VLCbtag;
     Float_t VLCdeltaeta;
     Float_t VLCdeltaphi;
     Float_t VLCdeltaR;
     Float_t Geneta;
     Float_t Genphi;
     Float_t Genpt;
     //Float_t Genbtag;
     Float_t Gendeltaeta;
     Float_t Gendeltaphi;
     Float_t GendeltaR;

     for(Long64_t entry=0; entry < nEntries; entry++){
	 tree_sig->GetEntry(entry);
	 tree_output->GetEntry(entry);
	 KTjet_size->GetBranch()->GetEntry(entry);
	 GenJet_size->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_size->GetBranch()->GetEntry(entry);
         
	 Int_t nKTjet = KTjet_size->GetValue();
	 Int_t nGenJet = GenJet_size->GetValue();
	 Int_t nVLCjet = VLCjetR05N4_size->GetValue();
	 
	 KTjet_eta->GetBranch()->GetEntry(entry);
         KTjet_phi->GetBranch()->GetEntry(entry);
	 KTjet_pt->GetBranch()->GetEntry(entry);
	 //KTjet_btag->GetBranch()->GetEntry(entry);
	 KTjet_deltaeta->GetBranch()->GetEntry(entry);
	 KTjet_deltaphi->GetBranch()->GetEntry(entry);

	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_pt->GetBranch()->GetEntry(entry);
	 //VLCjetR05N4_btag->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_deltaeta->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_deltaphi->GetBranch()->GetEntry(entry);

	 GenJet_eta->GetBranch()->GetEntry(entry);
	 GenJet_phi->GetBranch()->GetEntry(entry);
	 GenJet_pt->GetBranch()->GetEntry(entry);
	 //GenJet_btag->GetBranch()->GetEntry(entry);
	 GenJet_deltaeta->GetBranch()->GetEntry(entry);
	 GenJet_deltaphi->GetBranch()->GetEntry(entry);
         
	 //cout << "In the "<<entry+1<< "th event, among "<<nKTjet<<" KT, "<<nVLCjet <<" VLC, "<<nGenJet<<" Gen:" <<endl; 
	 
         Int_t nKTGen = 0;
	 Int_t nVLCGen = 0;
	 for(Long64_t ktentry=0; ktentry < nKTjet; ktentry++){
	     for(Long64_t genentry=0; genentry < nGenJet; genentry++){
		 KTeta = KTjet_eta->GetValue(ktentry);
         	 KTphi = KTjet_phi->GetValue(ktentry);
         	 KTpt = KTjet_pt->GetValue(ktentry);
         	 //KTbtag = KTjet_btag->GetValue(ktentry);
         	 KTdeltaeta = KTjet_deltaeta->GetValue(ktentry);
         	 KTdeltaphi = KTjet_deltaphi->GetValue(ktentry);
 	         //cout << "This is the "<< ktentry+1 << "th KTjet:"<<KTeta<<"|"<<KTphi<<endl; 	 
                 Geneta = GenJet_eta->GetValue(genentry);
         	 Genphi = GenJet_phi->GetValue(genentry);
         	 Genpt = GenJet_pt->GetValue(genentry);
         	 //Genbtag = GenJet_btag->GetValue(genentry);
         	 Gendeltaeta = GenJet_deltaeta->GetValue(genentry);
         	 Gendeltaphi = GenJet_deltaphi->GetValue(genentry);
 	         //cout << "This is the "<< genentry+1 << "th GenJet:"<<Geneta<<"|"<<Genphi<<endl;
		 Float_t DeltaR = TMath::Sqrt(pow((KTeta-Geneta),2)+pow((KTphi-Genphi),2));
		 if (DeltaR < 0.25){
		     //cout << ktentry+1 << "th  KT & "<<genentry+1 << "th Gen matched; ";
		     KTGeneta[nKTGen] = (KTeta - Geneta)/Geneta;
		     KTGenphi[nKTGen] = (KTphi - Genphi)/Genphi;
		     KTGenpt[nKTGen] = (KTpt - Genpt)/Genpt;
		     //KTGenbtag[nKTGen] = (KTbtag - Genbtag)/Genbtag;
		     Float_t KTdeltaR = TMath::Sqrt(pow((KTdeltaeta),2)+pow((KTdeltaphi),2));
		     Float_t GendeltaR = TMath::Sqrt(pow((Gendeltaeta),2)+pow((Gendeltaphi),2));
		     KTGendeltaR[nKTGen] = (KTdeltaR - GendeltaR)/GendeltaR;
		     nKTGen++;
		 }
	     }  
	 }
	 KTGenSize = nKTGen;
	 //cout <<nKTGen<< " KTjet & GenJet matched; ";
	 for(Long64_t vlcentry=0; vlcentry < nVLCjet; vlcentry++){
	     for(Long64_t genentry=0; genentry < nGenJet; genentry++){
		 VLCeta = VLCjetR05N4_eta->GetValue(vlcentry);
         	 VLCphi = VLCjetR05N4_phi->GetValue(vlcentry);
         	 VLCpt = VLCjetR05N4_pt->GetValue(vlcentry);
         	 //VLCbtag = VLCjetR05N4_btag->GetValue(vlcentry);
         	 VLCdeltaeta = VLCjetR05N4_deltaeta->GetValue(vlcentry);
         	 VLCdeltaphi = VLCjetR05N4_deltaphi->GetValue(vlcentry);
 	         //cout << "This is the "<< ktentry+1 << "th KTjet:"<<KTeta<<"|"<<KTphi<<endl; 	 
                 Geneta = GenJet_eta->GetValue(genentry);
         	 Genphi = GenJet_phi->GetValue(genentry);
         	 Genpt = GenJet_pt->GetValue(genentry);
         	 //Genbtag = GenJet_btag->GetValue(genentry);
         	 Gendeltaeta = GenJet_deltaeta->GetValue(genentry);
         	 Gendeltaphi = GenJet_deltaphi->GetValue(genentry);
 	         //cout << "This is the "<< genentry+1 << "th GenJet:"<<Geneta<<"|"<<Genphi<<endl;
		 Float_t DeltaR = TMath::Sqrt(pow((VLCeta-Geneta),2)+pow((VLCphi-Genphi),2));
		 if (DeltaR < 0.25){
		     //cout<< vlcentry+1 << "th VLC & "<<genentry+1 << "th Gen matched; ";
		     VLCR05N4Geneta[nVLCGen] = (VLCeta - Geneta)/Geneta;
		     VLCR05N4Genphi[nVLCGen] = (VLCphi - Genphi)/Genphi;
		     VLCR05N4Genpt[nVLCGen] = (VLCpt - Genpt)/Genpt;
		     //VLCR05N4Genbtag[nVLCGen] = (VLCbtag - Genbtag)/Genbtag;
                     Float_t VLCR05N4deltaR = TMath::Sqrt(pow((VLCdeltaeta),2)+pow((VLCdeltaphi),2));
		     Float_t GendeltaR = TMath::Sqrt(pow((Gendeltaeta),2)+pow((Gendeltaphi),2));
		     VLCR05N4GendeltaR[nVLCGen] = (VLCR05N4deltaR - GendeltaR)/GendeltaR;
		     nVLCGen++;
		 }	 
	     }  
	 }
	 //cout <<nVLCGen << " VLCjet & GenJet matched;" <<endl;
         VLCR05N4GenSize = nVLCGen;
	 tree_output->Fill();
     }
     tree_output->Write();
     TCanvas *mycanvas = new TCanvas("mycanvas","My Canvas",200,10,600,480);
     tree_output->Draw("KTGenpt");
     mycanvas->SaveAs("KTGenpt.png");
     tree_output->Draw("VLCR05N4Genpt");
     mycanvas->SaveAs("VLCGenpt.png");
     /*tree_output->Draw("KTGenbtag");
     mycanvas->SaveAs("KTGenbtag.png");
     tree_output->Draw("VLCR05N4Genbtag");
     mycanvas->SaveAs("VLCGenbtag.png");*/
     tree_output->Draw("KTGendeltaR");
     mycanvas->SaveAs("KTGenDeltaR.png");
     tree_output->Draw("VLCR05N4GendeltaR");
     mycanvas->SaveAs("VLCR05N4GenDeltaR.png");

     output->Close();
     file_sig->Close();
}

