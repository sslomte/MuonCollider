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
     TLeaf *AKTjet_size = tree_sig->GetLeaf("AKTjet_size");
     TLeaf *GenJet_size = tree_sig->GetLeaf("GenJet_size");
     TLeaf *VLCjetR05N4_size = tree_sig->GetLeaf("VLCjetR05N4_size");

     TLeaf *VLCjetR05N4_eta = tree_sig->GetLeaf("VLCjetR05N4.Eta");
     TLeaf *VLCjetR05N4_phi = tree_sig->GetLeaf("VLCjetR05N4.Phi");
     TLeaf *VLCjetR05N4_pt = tree_sig->GetLeaf("VLCjetR05N4.PT");
     TLeaf *VLCjetR05N4_deltaeta = tree_sig->GetLeaf("VLCjetR05N4.DeltaEta");
     TLeaf *VLCjetR05N4_deltaphi = tree_sig->GetLeaf("VLCjetR05N4.DeltaPhi");

     TLeaf *KTjet_eta = tree_sig->GetLeaf("KTjet.Eta");
     TLeaf *KTjet_phi = tree_sig->GetLeaf("KTjet.Phi");
     TLeaf *KTjet_pt= tree_sig->GetLeaf("KTjet.PT");
     TLeaf *KTjet_deltaeta = tree_sig->GetLeaf("KTjet.DeltaEta");
     TLeaf *KTjet_deltaphi = tree_sig->GetLeaf("KTjet.DeltaPhi");

     TLeaf *AKTjet_eta = tree_sig->GetLeaf("AKTjet.Eta");
     TLeaf *AKTjet_phi = tree_sig->GetLeaf("AKTjet.Phi");
     TLeaf *AKTjet_pt= tree_sig->GetLeaf("AKTjet.PT");
     TLeaf *AKTjet_deltaeta = tree_sig->GetLeaf("AKTjet.DeltaEta");
     TLeaf *AKTjet_deltaphi = tree_sig->GetLeaf("AKTjet.DeltaPhi");

     TLeaf *GenJet_eta = tree_sig->GetLeaf("GenJet.Eta");
     TLeaf *GenJet_phi = tree_sig->GetLeaf("GenJet.Phi");
     TLeaf *GenJet_pt = tree_sig->GetLeaf("GenJet.PT");
     TLeaf *GenJet_deltaeta = tree_sig->GetLeaf("GenJet.DeltaEta");
     TLeaf *GenJet_deltaphi = tree_sig->GetLeaf("GenJet.DeltaPhi");

     Int_t nEntries = tree_sig->GetEntries();

     TH1F *KTGenJetEta = new TH1F("KTGenJetEta", "KTGenJetEta", 100 , -2, 2); 
     TH1F *KTGenJetPhi = new TH1F("KTGenJetPhi", "KTGenJetPhi", 100 , -1, 1); 
     TH1F *KTGenJetPt = new TH1F("KTGenJetPt", "KTGenJetPt", 100 , -1, 1); 
     TH1F *KTGenJetDeltaR = new TH1F("KTGenJetDeltaR", "KTGenJetDeltaR", 100 , -6, 6); 
     TH1F *KTGenJetSize = new TH1F("KTGenJetSize", "KTGenJetSize", 8 , 0, 8); 

     TH1F *AKTGenJetEta = new TH1F("AKTGenJetEta", "AKTGenJetEta", 100 , -2, 2); 
     TH1F *AKTGenJetPhi = new TH1F("AKTGenJetPhi", "AKTGenJetPhi", 100 , -1, 1); 
     TH1F *AKTGenJetPt = new TH1F("AKTGenJetPt", "AKTGenJetPt", 100 , -1, 1); 
     TH1F *AKTGenJetDeltaR = new TH1F("AKTGenJetDeltaR", "AKTGenJetDeltaR", 100 , -6, 6); 
     TH1F *AKTGenJetSize = new TH1F("AKTGenJetSize", "AKTGenJetSize", 8 , 0, 8); 

     TH1F *VLCR05N4GenJetEta = new TH1F("VLCR05N4GenJetEta", "VLCR05N4GenJetEta", 100 , -2, 2); 
     TH1F *VLCR05N4GenJetPhi = new TH1F("VLCR05N4GenJetPhi", "VLCR05N4GenJetPhi", 100 , -1, 1); 
     TH1F *VLCR05N4GenJetPt = new TH1F("VLCR05N4GenJetPt", "VLCR05N4GenJetPt", 100 , -1, 1); 
     TH1F *VLCR05N4GenJetDeltaR = new TH1F("VLCR05N4GenJetDeltaR", "VLCR05N4GenJetDeltaR", 100 , -6, 6); 
     TH1F *VLCR05N4GenJetSize = new TH1F("VLCR05N4GenJetSize", "VLCR05N4GenJetSize", 8 , 0, 8); 


     Float_t KTGeneta;
     Float_t KTGenphi;
     Float_t KTGenpt;
     Float_t KTGendeltaR;
     Int_t KTGenSize;

     Float_t AKTGeneta;
     Float_t AKTGenphi;
     Float_t AKTGenpt;
     Float_t AKTGendeltaR;
     Int_t AKTGenSize;

     Float_t VLCR05N4Genphi;
     Float_t VLCR05N4Geneta;
     Float_t VLCR05N4Genpt;
     Float_t VLCR05N4GendeltaR;
     Int_t VLCR05N4GenSize;

     Float_t KTeta;
     Float_t KTphi;
     Float_t KTpt;
     Float_t KTdeltaeta;
     Float_t KTdeltaphi;
     Float_t KTdeltaR;

     Float_t AKTeta;
     Float_t AKTphi;
     Float_t AKTpt;
     Float_t AKTdeltaeta;
     Float_t AKTdeltaphi;
     Float_t AKTdeltaR;
          
     Float_t VLCeta;
     Float_t VLCphi;
     Float_t VLCpt;
     Float_t VLCdeltaeta;
     Float_t VLCdeltaphi;
     Float_t VLCdeltaR;
     
     Float_t Geneta;
     Float_t Genphi;
     Float_t Genpt;
     Float_t Gendeltaeta;
     Float_t Gendeltaphi;
     Float_t GendeltaR;

     for(Long64_t entry=0; entry < nEntries; entry++){
	 tree_sig->GetEntry(entry);
	 tree_output->GetEntry(entry);
	 KTjet_size->GetBranch()->GetEntry(entry);
	 AKTjet_size->GetBranch()->GetEntry(entry);
	 GenJet_size->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_size->GetBranch()->GetEntry(entry);
         
	 Int_t nKTjet = KTjet_size->GetValue();
	 Int_t nAKTjet = AKTjet_size->GetValue();
	 Int_t nGenJet = GenJet_size->GetValue();
	 Int_t nVLCjet = VLCjetR05N4_size->GetValue();
	 
	 KTjet_eta->GetBranch()->GetEntry(entry);
         KTjet_phi->GetBranch()->GetEntry(entry);
	 KTjet_pt->GetBranch()->GetEntry(entry);
	 KTjet_deltaeta->GetBranch()->GetEntry(entry);
	 KTjet_deltaphi->GetBranch()->GetEntry(entry);
 
	 AKTjet_eta->GetBranch()->GetEntry(entry);
         AKTjet_phi->GetBranch()->GetEntry(entry);
	 AKTjet_pt->GetBranch()->GetEntry(entry);
	 AKTjet_deltaeta->GetBranch()->GetEntry(entry);
	 AKTjet_deltaphi->GetBranch()->GetEntry(entry);

	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_pt->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_deltaeta->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_deltaphi->GetBranch()->GetEntry(entry);

	 GenJet_eta->GetBranch()->GetEntry(entry);
	 GenJet_phi->GetBranch()->GetEntry(entry);
	 GenJet_pt->GetBranch()->GetEntry(entry);
	 GenJet_deltaeta->GetBranch()->GetEntry(entry);
	 GenJet_deltaphi->GetBranch()->GetEntry(entry);
         
	 //cout << "In the "<<entry+1<< "th event, among "<<nKTjet<<" KT, "<<nVLCjet <<" VLC, "<<nGenJet<<" Gen:" <<endl; 
	 
         Int_t nKTGen = 0;
         Int_t nAKTGen = 0;
	 Int_t nVLCGen = 0;
	 for(Long64_t ktentry=0; ktentry < nKTjet; ktentry++){
	     for(Long64_t genentry=0; genentry < nGenJet; genentry++){
		 KTeta = KTjet_eta->GetValue(ktentry);
         	 KTphi = KTjet_phi->GetValue(ktentry);
         	 KTpt = KTjet_pt->GetValue(ktentry);
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
		     KTGeneta = (KTeta - Geneta)/Geneta;
		     KTGenphi = (KTphi - Genphi)/Genphi;
		     KTGenpt = (KTpt - Genpt)/Genpt;
		     Float_t KTdeltaR = TMath::Sqrt(pow((KTdeltaeta),2)+pow((KTdeltaphi),2));
		     Float_t GendeltaR = TMath::Sqrt(pow((Gendeltaeta),2)+pow((Gendeltaphi),2));
		     KTGendeltaR = (KTdeltaR - GendeltaR)/GendeltaR;
		     
		     KTGenJetEta->Fill(KTGeneta);
		     KTGenJetPhi->Fill(KTGenphi);
		     KTGenJetPt->Fill(KTGenpt);
		     KTGenJetDeltaR->Fill(KTGendeltaR);
		     nKTGen++;
		 }
	     }  
	 }
	 KTGenJetSize->Fill(nKTGen);

         for(Long64_t aktentry=0; aktentry < nAKTjet; aktentry++){
	     for(Long64_t genentry=0; genentry < nGenJet; genentry++){
		 AKTeta = AKTjet_eta->GetValue(aktentry);
         	 AKTphi = AKTjet_phi->GetValue(aktentry);
         	 AKTpt = AKTjet_pt->GetValue(aktentry);
         	 AKTdeltaeta = AKTjet_deltaeta->GetValue(aktentry);
         	 AKTdeltaphi = AKTjet_deltaphi->GetValue(aktentry);
 	         //cout << "This is the "<< ktentry+1 << "th KTjet:"<<KTeta<<"|"<<KTphi<<endl; 	 
                 Geneta = GenJet_eta->GetValue(genentry);
         	 Genphi = GenJet_phi->GetValue(genentry);
         	 Genpt = GenJet_pt->GetValue(genentry);
         	 //Genbtag = GenJet_btag->GetValue(genentry);
         	 Gendeltaeta = GenJet_deltaeta->GetValue(genentry);
         	 Gendeltaphi = GenJet_deltaphi->GetValue(genentry);
 	         //cout << "This is the "<< genentry+1 << "th GenJet:"<<Geneta<<"|"<<Genphi<<endl;
		 Float_t DeltaR = TMath::Sqrt(pow((AKTeta-Geneta),2)+pow((AKTphi-Genphi),2));
		 if (DeltaR < 0.25){
		     //cout << ktentry+1 << "th  KT & "<<genentry+1 << "th Gen matched; ";
		     AKTGeneta = (AKTeta - Geneta)/Geneta;
		     AKTGenphi = (AKTphi - Genphi)/Genphi;
		     AKTGenpt = (AKTpt - Genpt)/Genpt;
		     Float_t AKTdeltaR = TMath::Sqrt(pow((AKTdeltaeta),2)+pow((AKTdeltaphi),2));
		     Float_t GendeltaR = TMath::Sqrt(pow((Gendeltaeta),2)+pow((Gendeltaphi),2));
		     AKTGendeltaR = (AKTdeltaR - GendeltaR)/GendeltaR;
		     
		     AKTGenJetEta->Fill(AKTGeneta);
		     AKTGenJetPhi->Fill(AKTGenphi);
		     AKTGenJetPt->Fill(AKTGenpt);
		     AKTGenJetDeltaR->Fill(AKTGendeltaR);
		     nAKTGen++;
		 }
	     }  
	 }
	 AKTGenJetSize->Fill(nAKTGen);

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
		     VLCR05N4Geneta = (VLCeta - Geneta)/Geneta;
		     VLCR05N4Genphi = (VLCphi - Genphi)/Genphi;
		     VLCR05N4Genpt = (VLCpt - Genpt)/Genpt;
                     Float_t VLCR05N4deltaR = TMath::Sqrt(pow((VLCdeltaeta),2)+pow((VLCdeltaphi),2));
		     Float_t GendeltaR = TMath::Sqrt(pow((Gendeltaeta),2)+pow((Gendeltaphi),2));
		     VLCR05N4GendeltaR = (VLCR05N4deltaR - GendeltaR)/GendeltaR;
		     nVLCGen++;
		     VLCR05N4GenJetEta->Fill(VLCR05N4Geneta);
		     VLCR05N4GenJetPhi->Fill(VLCR05N4Genphi);
		     VLCR05N4GenJetPt->Fill(VLCR05N4Genpt);
		     VLCR05N4GenJetDeltaR->Fill(VLCR05N4GendeltaR);
		 }	 
	     }  
	 }
	 //cout <<nVLCGen << " VLCjet & GenJet matched;" <<endl;
         VLCR05N4GenJetSize->Fill(nVLCGen);
     }
     KTGenJetEta->Write();
     KTGenJetPhi->Write();
     KTGenJetPt->Write();
     KTGenJetDeltaR->Write();
     KTGenJetSize->Write();

     AKTGenJetEta->Write();
     AKTGenJetPhi->Write();
     AKTGenJetPt->Write();
     AKTGenJetDeltaR->Write();
     AKTGenJetSize->Write();

     VLCR05N4GenJetEta->Write();
     VLCR05N4GenJetPhi->Write();
     VLCR05N4GenJetPt->Write();
     VLCR05N4GenJetDeltaR->Write();
     VLCR05N4GenJetSize->Write();

     tree_output->Write();
     
     TCanvas *mycanvas = new TCanvas("mycanvas","My Canvas",200,10,600,480);
     KTGenJetPt->Draw();
     mycanvas->SaveAs("KTGenpt.png");
     VLCR05N4GenJetPt->Draw();
     mycanvas->SaveAs("VLCGenpt.png");
     KTGenJetDeltaR->Draw();
     mycanvas->SaveAs("KTGenDeltaR.png");
     VLCR05N4GenJetDeltaR->Draw();
     mycanvas->SaveAs("VLCR05N4GenDeltaR.png");
     KTGenJetSize->Draw();
     mycanvas->SaveAs("KTGenSize.png");
     VLCR05N4GenJetSize->Draw();
     mycanvas->SaveAs("VLCR05N4GenSize.png");

     AKTGenJetPt->Draw();
     mycanvas->SaveAs("AKTGenpt.png");
     AKTGenJetDeltaR->Draw();
     mycanvas->SaveAs("AKTGenDeltaR.png");
     AKTGenJetSize->Draw();
     mycanvas->SaveAs("AKTGenSize.png");

     output->Close();
     file_sig->Close();
}

