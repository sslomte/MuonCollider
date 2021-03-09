import ROOT
from ROOT import TFile, gROOT, TCanvas, gSystem
import numpy as np
gSystem.Load("../../libDelphes.so")

print "...import signal and background file"
input_sig = "../../delphes_dhiggs_bkgvvbbz.root"
file_sig = ROOT.TFile(input_sig)
tree_sig = ROOT.gROOT.FindObject("Delphes")
myCanvas = TCanvas()
branch_list = ["KTjet","Particle", "MissingET", "Track"]
vars_toprint = ["Phi", "PT", "Eta", "E", "Mass", "Charge"]
keys = tree_sig.GetListOfLeaves()
branchlist = tree_sig.GetListOfLeaves()
for branch in branch_list:
	for var in vars_toprint:
		name = branch+"."+var
		print "print "+name+"..."
		if (name in keys):
			tree_sig.Draw(name)
			myCanvas.Draw()
			myCanvas.SaveAs(branch+"_"+var+".png")
		else: 
			print var+" is not in "+branch
