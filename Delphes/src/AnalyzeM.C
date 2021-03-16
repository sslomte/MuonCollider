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
     //TChain chain("Delphes");
     //chain.Add(inputFile);
     //ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
     TFile *file_sig = new TFile(inputFile);
     TFile *output = new TFile(outputFile, "recreate");
     TTree *tree_sig = (TTree*)file_sig->Get("Delphes");
     TTree *tree_output = new TTree("tree_output","Delphes");

     TLeaf *VLCjetR05N4_size = tree_sig->GetLeaf("VLCjetR05N4_size");
     TLeaf *VLCjetR05N4_eta = tree_sig->GetLeaf("VLCjetR05N4.Eta");
     TLeaf *VLCjetR05N4_phi = tree_sig->GetLeaf("VLCjetR05N4.Phi");
     TLeaf *VLCjetR05N4_pt = tree_sig->GetLeaf("VLCjetR05N4.PT");
     TLeaf *VLCjetR05N4_mass = tree_sig->GetLeaf("VLCjetR05N4.Mass");

     TLeaf *GenJet_size = tree_sig->GetLeaf("GenJet_size");
     TLeaf *GenJet_eta = tree_sig->GetLeaf("GenJet.Eta");
     TLeaf *GenJet_phi = tree_sig->GetLeaf("GenJet.Phi");

     Int_t nEntries = tree_sig->GetEntries();

     TH1D *VLCR05N4Mass1 = new TH1D("VLCR05N4Mass1", "VLCR05N4Mass1", 100 , 0, 500); 
     TH1D *VLCR05N4Mass2 = new TH1D("VLCR05N4Mass2", "VLCR05N4Mass2", 100 , 0, 500); 

     Double_t VLC1eta;
     Double_t VLC1phi;
     Double_t VLC1pt;
     Double_t VLC1mass;
     Double_t VLC2eta;
     Double_t VLC2phi;
     Double_t VLC2pt;
     Double_t VLC2mass;
     Double_t VLCR05N4pairmass;
     Double_t Gen1eta;
     Double_t Gen1phi;
     Double_t Gen2eta;
     Double_t Gen2phi;

     Int_t pair1jet1entry;
     Int_t pair1jet2entry;
     Int_t pair2jet1entry;
     Int_t pair2jet2entry;

     cout << "start";
     for(Long64_t entry=0; entry < nEntries; entry++){
	 tree_sig->GetEntry(entry);
	 tree_output->GetEntry(entry);
	 VLCjetR05N4_size->GetBranch()->GetEntry(entry);
         
	 Int_t nVLCjet = VLCjetR05N4_size->GetValue();
	 Int_t nGenJet = GenJet_size->GetValue();

	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_phi->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_pt->GetBranch()->GetEntry(entry);
	 VLCjetR05N4_mass->GetBranch()->GetEntry(entry);
	 Double_t VLCR05N4pair1Mass = 0;
	 Double_t VLCR05N4subpairmass = 0;

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
		 VLCR05N4pairmass = 0;

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

	 bool VLC1flag = false;
	 bool VLC2flag = false;
         Float_t jet1DeltaR = 100;
         Float_t jet2DeltaR = 100;
	 Int_t jet1entry;

	 for(Int_t gen1entry=0; gen1entry < nGenJet; gen1entry++){
  	     Gen1eta = GenJet_eta->GetValue(gen1entry);
	     Gen1phi = GenJet_phi->GetValue(gen1entry);
	     Float_t jet1DeltaRtmp = TMath::Sqrt(pow((VLC1eta-Gen1eta),2)+pow((VLC1phi-Gen1phi),2));
	     //cout << jet1DeltaRtmp << "|";
	     if(jet1DeltaRtmp < jet1DeltaR){
	         jet1DeltaR = jet1DeltaRtmp;
		 jet1entry = gen1entry;
		 if (jet1DeltaR < 0.25){
		     VLC1flag = true;
	             //cout << "the first jet in the "<<entry<< "th event pass the check."<<endl;	 
		 }
	     }  
         }
         for(Int_t gen2entry=0; gen2entry < nGenJet; gen2entry++){
  	     Gen2eta = GenJet_eta->GetValue(gen2entry);
	     Gen2phi = GenJet_phi->GetValue(gen2entry);
	     Float_t jet2DeltaRtmp = TMath::Sqrt(pow((VLC2eta-Gen2eta),2)+pow((VLC2phi-Gen2phi),2));
	     if (gen2entry != jet1entry){
                 if(jet2DeltaRtmp < jet2DeltaR){
	             jet2DeltaR = jet2DeltaRtmp;
		     Int_t jet2entry = gen2entry;
                     if (jet2DeltaR < 0.25){
		         VLC2flag = true;
	                 //cout << "the second jet in the "<<entry<< "th event pass the check."<<endl;	 
			   
		     }
	         }
	     }
         }
	 
         Double_t VLCR05N4pair2Mass = 0;
         TLorentzVector h2;
         TLorentzVector jet1;
         TLorentzVector jet2;
 	 jet1.SetPtEtaPhiM(VLC1pt, VLC1eta, VLC1phi,VLC1mass);
	 jet2.SetPtEtaPhiM(VLC2pt, VLC2eta, VLC2phi,VLC1mass);
	 h2=jet1+jet2;
	 VLCR05N4pair2Mass = h2.Mag();
         //VLCR05N4pair2Mass = TMath::Sqrt(2*VLC1pt*VLC2pt*(TMath::CosH(VLC1eta-VLC2eta)-TMath::Cos(VLC1phi-VLC2phi)));
     //VLCR05N4Mass1->Fill(VLCR05N4pair1Mass);
     //VLCR05N4Mass2->Fill(VLCR05N4pair2Mass);
	 
         if(VLC1flag==true and VLC2flag==true){
	     VLCR05N4Mass1->Fill(VLCR05N4pair1Mass);
             VLCR05N4Mass2->Fill(VLCR05N4pair2Mass);
	 } else {
	     cout << "the second jet pair in the "<<entry<< "th event failed check."<<endl;	 
         }
     }
     TF1 *jetpair1gausfit = new TF1("jetpair1gausfit", "gaus",20,220);
     TF1 *jetpair2gausfit = new TF1("jetpair2gausfit", "gaus",0,140);
     
     TCanvas *mycanvas = new TCanvas("mycanvas","My Canvas",200,10,600,480);
     VLCR05N4Mass1->Fit("jetpair1gausfit","R");
     mycanvas->SaveAs("VLCR05N4pair1Mass.png");
     VLCR05N4Mass2->Fit("jetpair2gausfit","R");
     mycanvas->SaveAs("VLCR05N4pair2Mass.png");
     VLCR05N4Mass1->Write();
     VLCR05N4Mass2->Write();

     tree_output->Write();


     output->Close();
     file_sig->Close();
}

