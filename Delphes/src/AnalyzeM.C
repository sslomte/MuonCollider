//usage: root -l AnalyzeM.C\(\"inputfile.root\"\,\"outputfile.root\"\)
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
#include <TLorentzVector.h>

void AnalyzeM(const char *inputFile, const char *outputFile){
     gSystem->Load("libDelphes.so");
     TChain chain("Delphes");
     chain.Add(inputFile);
     ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
     TFile *file_sig = new TFile(inputFile);
     TFile *output = new TFile(outputFile, "recreate");
     TTree *tree_sig = (TTree*)file_sig->Get("Delphes");
     TTree *tree_output = new TTree("tree_output","Delphes");
     TLeaf *VLCjetR05N4_size = tree_sig->GetLeaf("VLCjetR05N4_size");

     TLeaf *VLCjetR05N4_eta = tree_sig->GetLeaf("VLCjetR05N4.Eta");
     TLeaf *VLCjetR05N4_phi = tree_sig->GetLeaf("VLCjetR05N4.Phi");
     TLeaf *VLCjetR05N4_pt = tree_sig->GetLeaf("VLCjetR05N4.PT");
     TLeaf *VLCjetR05N4_mass = tree_sig->GetLeaf("VLCjetR05N4.Mass");

     Int_t nEntries = tree_sig->GetEntries();

     TH1D *VLCR05N4Mass1 = new TH1D("VLCR05N4Mass1", "VLCR05N4Mass1", 100 , 0, 1000); 
     TH1D *VLCR05N4Mass2 = new TH1D("VLCR05N4Mass2", "VLCR05N4Mass2", 100 , 0, 1000); 

     Double_t VLC1eta;
     Double_t VLC1phi;
     Double_t VLC1pt;
     Double_t VLC1mass;
     Double_t VLC2eta;
     Double_t VLC2phi;
     Double_t VLC2pt;
     Double_t VLC2mass;
     Double_t VLCR05N4pairmass;
     Double_t VLCR05N4pair1Mass = 0;
     Double_t VLCR05N4pair2Mass = 0;
     Int_t pair1jet1entry;
     Int_t pair1jet2entry;
     Int_t pair2jet1entry;
     Int_t pair2jet2entry;

     for(Long64_t entry=0; entry < nEntries; entry++){
	 tree_sig->GetEntry(entry);
	 tree_output->GetEntry(entry);
	 VLCjetR05N4_size->GetBranch()->GetEntry(entry);
         
	 Int_t nVLCjet = VLCjetR05N4_size->GetValue();

	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_pt->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_mass->GetBranch()->GetEntry(entry);
	 
	 for(Int_t vlc1entry=0; vlc1entry < nVLCjet; vlc1entry++){
	     for(Int_t vlc2entry=vlc1entry+1; vlc2entry < nVLCjet; vlc2entry++){
		 VLC1eta = VLCjetR05N4_eta->GetValue(vlc1entry);
         	 VLC1phi = VLCjetR05N4_phi->GetValue(vlc1entry);
         	 VLC1pt = VLCjetR05N4_pt->GetValue(vlc1entry);
         	 VLC1mass = VLCjetR05N4_mass->GetValue(vlc1entry);
                 VLC2eta = VLCjetR05N4_eta->GetValue(vlc2entry);
         	 VLC2phi = VLCjetR05N4_phi->GetValue(vlc2entry);
         	 VLC2pt = VLCjetR05N4_pt->GetValue(vlc2entry);
         	 VLC2mass = VLCjetR05N4_mass->GetValue(vlc2entry);

                 TLorentzVector h1;
                 TLorentzVector jet1;
                 TLorentzVector jet2;
		 jet1.SetPtEtaPhiM(VLC1pt, VLC1eta, VLC1phi,VLC1mass);
		 jet2.SetPtEtaPhiM(VLC2pt, VLC2eta, VLC2phi,VLC1mass);
		 h1=jet1+jet2;
		 VLCR05N4pairmass = h1.Mag();
                 //VLCR05N4pairmass = TMath::Sqrt(2*VLC1pt*VLC2pt*(TMath::CosH(VLC1eta-VLC2eta)-TMath::Cos(VLC1phi-VLC2phi)));
		 if(abs(125 - VLCR05N4pairmass) < abs(125 -VLCR05N4pair1Mass)){
		     VLCR05N4pair1Mass = VLCR05N4pairmass;
		     pair1jet1entry = vlc1entry;
		     pair1jet2entry = vlc2entry;
		 }		 	 
	     }  
	 }

	 for(Int_t vlcentry=0; vlcentry < nVLCjet; vlcentry++){
	     if((vlcentry!=pair1jet1entry) and (vlcentry!=pair1jet2entry)){
                 pair2jet1entry = vlcentry;
		 break;
	     }
	 }
         for(Int_t vlcentry=pair2jet1entry+1; vlcentry < nVLCjet; vlcentry++){
	     if(((vlcentry!=pair1jet1entry) and (vlcentry!=pair1jet2entry)) and (vlcentry!=pair2jet1entry)){
	         pair2jet2entry = vlcentry;	     
	     }
	 }
         VLC1eta = VLCjetR05N4_eta->GetValue(pair2jet1entry);
         VLC1phi = VLCjetR05N4_phi->GetValue(pair2jet1entry);
         VLC1pt = VLCjetR05N4_pt->GetValue(pair2jet1entry);
         VLC1mass = VLCjetR05N4_mass->GetValue(pair2jet1entry);
         VLC2eta = VLCjetR05N4_eta->GetValue(pair2jet2entry);
         VLC2phi = VLCjetR05N4_phi->GetValue(pair2jet2entry);
         VLC2pt = VLCjetR05N4_pt->GetValue(pair2jet2entry);
         VLC2mass = VLCjetR05N4_mass->GetValue(pair2jet2entry);

         TLorentzVector h2;
         TLorentzVector jet1;
         TLorentzVector jet2;
 	 jet1.SetPtEtaPhiM(VLC1pt, VLC1eta, VLC1phi,VLC1mass);
	 jet2.SetPtEtaPhiM(VLC2pt, VLC2eta, VLC2phi,VLC1mass);
	 h2=jet1+jet2;
	 VLCR05N4pair2Mass = h2.Mag();
         //VLCR05N4pair2Mass = TMath::Sqrt(2*VLC1pt*VLC2pt*(TMath::CosH(VLC1eta-VLC2eta)-TMath::Cos(VLC1phi-VLC2phi)));
	 
	 VLCR05N4Mass1->Fill(VLCR05N4pair1Mass);
         VLCR05N4Mass2->Fill(VLCR05N4pair2Mass);

	 //cout << pair1jet1entry <<"|"<<pair1jet2entry<<"|"<<pair2jet1entry<<"|"<<pair2jet2entry<<endl;
     }
     VLCR05N4Mass1->Write();
     VLCR05N4Mass2->Write();

     tree_output->Write();
     TCanvas *mycanvas = new TCanvas("mycanvas","My Canvas",200,10,600,480);
     VLCR05N4Mass1->Draw();
     mycanvas->SaveAs("VLCR05N4pair1Mass.png");
     VLCR05N4Mass2->Draw();
     mycanvas->SaveAs("VLCR05N4pair2Mass.png");
     
     output->Close();
     file_sig->Close();
}

