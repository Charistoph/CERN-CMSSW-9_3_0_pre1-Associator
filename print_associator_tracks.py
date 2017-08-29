import os
import sys
import ROOT
from math import sqrt,pi

# open root file & tree
filename = "output_gsf_associator.root"
tf = ROOT.TFile(filename)
#tf.ls()
tree_dir = tf.Get("MyTrackAssociator")
tree = tree_dir.Get("track_associator_tree")
#tree_dir.ls()

# create output folder
outputpath = filename.replace(".root", "")
#namestr = filename.replace(".root", "")
#outputpath = os.path.join('output_' + namestr)
if (os.path.exists(outputpath)==False):
    os.makedirs(outputpath)
    print "path created", outputpath

print "path and file routine complete\n"

# define histograms
assoc_track_pt = ROOT.TH1F("assoc_track_pt", "assoc_track_pt", 100, 0.0, 100.0)
assoc_track_phi = ROOT.TH1F("assoc_track_phi", "assoc_track_phi", 3, -3, 3)
assoc_track_eta = ROOT.TH1F("assoc_track_eta", "assoc_track_eta", 3, -3, 3)
assoc_track_numberOfValidHits = ROOT.TH1F("assoc_track_numberOfValidHits", "assoc_track_numberOfValidHits", 100, 0.0, 100.0)
all_track_quality = ROOT.TH1F("all_track_quality", "all_track_quality", 2, -2.0, 2.0)
all_track_pt = ROOT.TH1F("all_track_pt", "all_track_pt", 100, 0.0, 100.0)
all_track_phi = ROOT.TH1F("all_track_phi", "all_track_phi", 3, -3, 3)
all_track_eta = ROOT.TH1F("all_track_eta", "all_track_eta", 3, -3, 3)
all_track_numberOfValidHits = ROOT.TH1F("all_track_eta", "all_track_numberOfValidHits", 100, 0.0, 100.0)

# loop over all events in tree
for event in tree:
    if event.assoc_track[7] == 1:
        assoc_track_pt.Fill(event.assoc_track[0])
        assoc_track_phi.Fill(event.assoc_track[1])
        assoc_track_eta.Fill(event.assoc_track[2])
        assoc_track_numberOfValidHits.Fill(event.assoc_track[6])
    all_track_quality.Fill(event.assoc_track[7])
    all_track_pt.Fill(event.assoc_track[0])
    all_track_phi.Fill(event.assoc_track[1])
    all_track_eta.Fill(event.assoc_track[2])
    all_track_numberOfValidHits.Fill(event.assoc_track[6])

c = ROOT.TCanvas( "c", "c", 800, 800)

# save histos
assoc_track_pt.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/assoc_track_pt.png")

assoc_track_phi.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/assoc_track_phi.png")

assoc_track_eta.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/assoc_track_eta.png")

assoc_track_numberOfValidHits.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/assoc_track_numberOfValidHits.png")

all_track_quality.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/all_track_quality.png")

all_track_pt.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/all_track_pt.png")

all_track_phi.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/all_track_phi.png")

all_track_eta.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/all_track_eta.png")

all_track_numberOfValidHits.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/all_track_numberOfValidHits.png")

print "Histos printed.\n"

