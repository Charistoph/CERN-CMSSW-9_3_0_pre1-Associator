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

# loop over all events in tree
for event in tree:
    gsf_track_pt.Fill(event.gsf_track[0])

c = ROOT.TCanvas( "c", "c", 800, 800)

# save histos
gsf_track_pt.Draw()
ROOT.gPad.Update()
c.SaveAs(outputpath + "/gsf_track_pt.png")

