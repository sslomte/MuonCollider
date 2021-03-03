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

void JetsMerge(const char *inputFile){
     gSystem->Load("libDelphes.so");
     TChain chain("Delphes");
     chain.Add(inputFile);
     ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
     TFile *file_sig = new TFile(inputFile);
     TFile *output = new TFile("delphes_dhiggs_sig_matchjets.root", "recreate");
     TTree *tree_sig = (TTree*)file_sig->Get("Delphes");
     TTree *tree_output = new TTree("tree_output","Delphes");
     TBranch *KTjet = tree_sig->GetBranch("KTjet");
     TBranch *GenJet = tree_sig->GetBranch("GenJet");
     TBranch *VLCjetR05N4 = tree_sig->GetBranch("VLCjetR05N4");
     TLeaf *KTjet_size = tree_sig->GetLeaf("KTjet_size");
     TLeaf *GenJet_size = tree_sig->GetLeaf("GenJet_size");
     TLeaf *VLCjetR05N4_size = tree_sig->GetLeaf("VLCjetR05N4_size");
     TLeaf *VLCjetR05N4_eta = tree_sig->GetLeaf("VLCjetR05N4.Eta");
     TLeaf *VLCjetR05N4_phi = tree_sig->GetLeaf("VLCjetR05N4.Phi");
     TLeaf *KTjet_eta = tree_sig->GetLeaf("KTjet.Eta");
     TLeaf *KTjet_phi = tree_sig->GetLeaf("KTjet.Phi");
     TLeaf *GenJet_eta = tree_sig->GetLeaf("GenJet.Eta");
     TLeaf *GenJet_phi = tree_sig->GetLeaf("GenJet.Phi");

     Int_t nEntries = tree_sig->GetEntries();
     
     Float_t KTGeneta[6];
     Float_t KTGenphi[6];
     Int_t KTGenSize;
     
     Float_t VLCR05N4Genphi[4];
     Float_t VLCR05N4Geneta[4];
     Int_t VLCR05N4GenSize;

     tree_output->Branch("KTGenJetEta", &KTGeneta, "KTGeneta/F");
     tree_output->Branch("KTGenJetPhi", &KTGenphi, "KTGenphi/F");
     tree_output->Branch("KTGenJetSize", &KTGenSize, "KTGenSize/I");
     
     tree_output->Branch("VLCR05N4GenJetEta", &VLCR05N4Geneta, "VLCR05N4Geneta/F");
     tree_output->Branch("VLCR05N4GenJetPhi", &VLCR05N4Genphi, "VLCR05N4Genphi/F");
     tree_output->Branch("VLCR05N4GenJetSize", &VLCR05N4GenSize, "VLCR05N4GenSize/I");

     Float_t KTeta;
     Float_t KTphi;
     Float_t VLCeta;
     Float_t VLCphi;
     Float_t Geneta;
     Float_t Genphi;

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
	 
	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 
	 GenJet_eta->GetBranch()->GetEntry(entry);
	 GenJet_phi->GetBranch()->GetEntry(entry);
         
	 //cout << "In the "<<entry+1<< "th event, among "<<nKTjet<<" KT, "<<nVLCjet <<" VLC, "<<nGenJet<<" Gen:" <<endl; 
	 
         Int_t nKTGen = 0;
	 Int_t nVLCGen = 0;
	 for(Long64_t ktentry=0; ktentry < nKTjet; ktentry++){
	     for(Long64_t genentry=0; genentry < nGenJet; genentry++){
		 KTeta = KTjet_eta->GetValue(ktentry);
         	 KTphi = KTjet_phi->GetValue(ktentry);
 	         //cout << "This is the "<< ktentry+1 << "th KTjet:"<<KTeta<<"|"<<KTphi<<endl; 	 
                 Geneta = GenJet_eta->GetValue(genentry);
         	 Genphi = GenJet_phi->GetValue(genentry);
 	         //cout << "This is the "<< genentry+1 << "th GenJet:"<<Geneta<<"|"<<Genphi<<endl;
		 Float_t DeltaR = TMath::Sqrt(pow((KTeta-Geneta),2)+pow((KTphi-Genphi),2));
		 if (DeltaR < 0.25){
		     //cout << ktentry+1 << "th  KT & "<<genentry+1 << "th Gen matched; ";
		     KTGeneta[nKTGen] = (KTeta + Geneta)/2;
		     KTGenphi[nKTGen] = (KTphi + Genphi)/2;
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
 	         //cout << "This is the "<< ktentry+1 << "th KTjet:"<<KTeta<<"|"<<KTphi<<endl; 	 
                 Geneta = GenJet_eta->GetValue(genentry);
         	 Genphi = GenJet_phi->GetValue(genentry);
 	         //cout << "This is the "<< genentry+1 << "th GenJet:"<<Geneta<<"|"<<Genphi<<endl;
		 Float_t DeltaR = TMath::Sqrt(pow((VLCeta-Geneta),2)+pow((VLCphi-Genphi),2));
		 if (DeltaR < 0.25){
		     //cout<< vlcentry+1 << "th VLC & "<<genentry+1 << "th Gen matched; ";
		     VLCR05N4Geneta[nVLCGen] = (VLCeta + Geneta)/2;
		     VLCR05N4Genphi[nVLCGen] = (VLCphi + Genphi)/2;
		     nVLCGen++;
		 }	 
	     }  
	 }
	 //cout <<nVLCGen << " VLCjet & GenJet matched;" <<endl;
         VLCR05N4GenSize = nVLCGen;
	 tree_output->Fill();
     }
     tree_output->Write();
     output->Close();
     file_sig->Close();
}

