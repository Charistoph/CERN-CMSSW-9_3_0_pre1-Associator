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
namestr = filename.replace(".root", "")
outputpath = os.path.join('output_' + namestr)
if (os.path.exists(outputpath)==False):
    os.makedirs(outputpath)
    print "path created", outputpath

print "path and file routine complete\n"

# define histograms
gsf_track_pt   = ROOT.TH1F("gsf_track_pt", "gsf_track_pt", 100, 0.0, 100.0)
gsf_track_phi   = ROOT.TH1F("gsf_track_phi", "gsf_track_phi", 100, 0.0, 100.0)
gsf_track_eta   = ROOT.TH1F("gsf_track_eta", "gsf_track_eta", 100, 0.0, 100.0)
sts_track_pt   = ROOT.TH1F("sts_track_pt", "sts_track_pt", 100, 0.0, 100.0)
sts_track_phi   = ROOT.TH1F("sts_track_phi", "sts_track_phi", 100, 0.0, 100.0)
sts_track_eta   = ROOT.TH1F("sts_track_eta", "sts_track_eta", 100, 0.0, 100.0)

# loop over all events in tree
for event in tree:
    gsf_track_pt.Fill(event.gsf_track[0])
    gsf_track_phi.Fill(event.gsf_track[1])
    gsf_track_eta.Fill(event.gsf_track[2])
    sts_track_pt.Fill(event.sts_track[0])
    sts_track_phi.Fill(event.sts_track[1])
    sts_track_eta.Fill(event.sts_track[2])

c = ROOT.TCanvas( "c", "c", 800, 800)

# save histos
gsf_track_pt.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/gsf_track_pt.png")

gsf_track_phi.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/gsf_track_phi.png")

gsf_track_eta.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/gsf_track_eta.png")

sts_track_pt.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/sts_track_pt.png")

sts_track_phi.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/sts_track_phi.png")

sts_track_eta.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/sts_track_eta.png")

print "Histos printed.\n"

