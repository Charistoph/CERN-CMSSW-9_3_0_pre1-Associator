import os
import sys
import ROOT
from math import sqrt,pi

#-------------------------------------------------------------------------------
# Define Funcitons

# define Histo Draw
def HistoDraw(track_name,track,variable_name):
    track.GetMinimum(0.)
    track.Draw()
    ROOT.gPad.Update()
    c.SaveAs(outputpath + "/" + track_name + "_" + variable_name + ".png")

# define Histo Fill
def HistoFill(event,pt,phi,eta,numberOfValidHits):
    pt.Fill(event.assoc_track[0])
    phi.Fill(event.assoc_track[1])
    eta.Fill(event.assoc_track[2])
    numberOfValidHits.Fill(event.assoc_track[6])

# define histograms
# assoc_para:
#    0 = seed and track not associated
#    1 = seed associated, track not associated
#    2 = seed and track associated
def MakeHisto(track_name,assoc_para):

# n pins, low, up
    pt = ROOT.TH1F( track_name + "_pt", track_name + "_pt", 100, 0.0, 100.0)
    phi = ROOT.TH1F( track_name + "_phi", track_name + "_phi", 25, -3.2, 3.2)
    eta = ROOT.TH1F( track_name + "_eta", track_name + "_eta", 25, -2.5, 2.5)
    numberOfValidHits = ROOT.TH1F( track_name + "_numberOfValidHits", track_name + "_numberOfValidHits", 25, 0.0, 25.0)

    print "Histos initiated."

    if assoc_para == 0:
        print "Making histogramms for non associated tracks"
        for event in tree:
            if (event.assoc_track[7] == -1 and event.assoc_track[8] == -1):
# Calling Fill Function
                HistoFill(event,pt,phi,eta,numberOfValidHits)

    if assoc_para == 1:
        print "Making histogramms for seed associated, track non associated tracks"
        for event in tree:
            if (event.assoc_track[7] == 1 and event.assoc_track[8] == -1):
# Calling Fill Function
                HistoFill(event,pt,phi,eta,numberOfValidHits)

    if assoc_para == 2:
        print "Making histogramms for seed and track associated tracks"
        for event in tree:
            if (event.assoc_track[7] == 1 and event.assoc_track[8] == 1):
# Calling Fill Function
                HistoFill(event,pt,phi,eta,numberOfValidHits)

    print "Histos filled."

# Calling Draw Function
    HistoDraw(track_name,pt,"pt")
    HistoDraw(track_name,phi,"phi")
    HistoDraw(track_name,eta,"eta")
    HistoDraw(track_name,numberOfValidHits,"numberOfValidHits")

    print "Histos saved.\n"

#-------------------------------------------------------------------------------
# main code

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
    print "\npath created", outputpath

print "\nPath and file routine complete.\n"

# Prepare Histo Canvas
c = ROOT.TCanvas( "c", "c", 800, 800)

# call Make Histo function
MakeHisto("non_assoc_track",0)
#MakeHisto("seed_assoc_track",1)
#MakeHisto("all_assoc_track",2)

print "All Histos printed."

