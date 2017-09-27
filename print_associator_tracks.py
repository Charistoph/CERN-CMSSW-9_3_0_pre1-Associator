import os
import sys
import ROOT
from math import sqrt,pi

#-------------------------------------------------------------------------------
# Define Funcitons

# Draw One Histo
def HistoDraw(track_name,track,variable_name):
    track.GetMinimum(0.)
    track.Draw()
    ROOT.gPad.Update()
    c.SaveAs(outputpath + "/" + track_name + "_" + variable_name + ".png")

def FillStandardHistos(branch_variable,event,pt,phi,eta,numberOfValidHits):
    pt.Fill(eval("event." + branch_variable + "[0]"))
    phi.Fill(eval("event." + branch_variable + "[1]"))
    eta.Fill(eval("event." + branch_variable + "[2]"))
    numberOfValidHits.Fill(eval("event." + branch_variable + "[6]"))

def DrawStandardHistos(track_name,pt,phi,eta,numberOfValidHits):
    HistoDraw(track_name,pt,"pt")
    HistoDraw(track_name,phi,"phi")
    HistoDraw(track_name,eta,"eta")
    HistoDraw(track_name,numberOfValidHits,"numberOfValidHits")

# define assoc para histograms
# assoc_para:
#    0 = seed and track not associated
#    1 = seed associated, track not associated
#    2 = seed and track associated
def MakeAssocParaHistos(branch_variable,track_name,assoc_para):

# parameters: n pins, low, up
    pt = ROOT.TH1F( track_name + "_pt", track_name + "_pt", 100, 0.0, 100.0)
    phi = ROOT.TH1F( track_name + "_phi", track_name + "_phi", 25, -3.2, 3.2)
    eta = ROOT.TH1F( track_name + "_eta", track_name + "_eta", 25, -2.5, 2.5)
    numberOfValidHits = ROOT.TH1F( track_name + "_numberOfValidHits", track_name + "_numberOfValidHits", 25, 0.0, 25.0)

    print "AssocPara Histos initiated."

    if assoc_para == 0:
        print "Making histogramms for non associated tracks"
        for event in tree:
            if (eval("event." + branch_variable + "[7]") == -1 and eval("event." + branch_variable + "[8]") == -1):
                FillStandardHistos(branch_variable,event,pt,phi,eta,numberOfValidHits)

    if assoc_para == 1:
        print "Making histogramms for seed associated, track non associated tracks"
        for event in tree:
            if (eval("event." + branch_variable + "[7]") == 1 and eval("event." + branch_variable + "[8]") == -1):
                FillStandardHistos(branch_variable,event,pt,phi,eta,numberOfValidHits)

    if assoc_para == 2:
        print "Making histogramms for seed and track associated tracks"
        for event in tree:
            if (eval("event." + branch_variable + "[7]") == 1 and eval("event." + branch_variable + "[8]") == 1):
                FillStandardHistos(branch_variable,event,pt,phi,eta,numberOfValidHits)

    print "Histos filled."
    DrawStandardHistos(track_name,pt,phi,eta,numberOfValidHits)
    print "Histos saved.\n"

def MakeQualHistos(branch_variable,track_name):

# n pins, low, up
    seed_qual = ROOT.TH1F( track_name + "_seed_qual", track_name + "_seed_qual", 3, -1.0, 2.0)
    track_qual = ROOT.TH1F( track_name + "_track_qual", track_name + "_track_qual", 3, -1.0, 2.0)

    print "Qual Histos initiated."

    for event in tree:
        seed_qual.Fill(eval("event." + branch_variable + "[7]"))
        track_qual.Fill(eval("event." + branch_variable + "[8]"))

    print "Histos filled."

# Calling Draw Function
    HistoDraw(track_name,seed_qual,"seed_qual")
    HistoDraw(track_name,track_qual,"track_qual")

    print "Histos saved.\n"

def MakeHistoIfNonZero(branch_variable,track_name):

    print "NonZero Histo initiated."

# parameters: n pins, low, up
    pt = ROOT.TH1F( track_name + "_pt", track_name + "_pt", 100, 0.0, 100.0)
    phi = ROOT.TH1F( track_name + "_phi", track_name + "_phi", 25, -3.2, 3.2)
    eta = ROOT.TH1F( track_name + "_eta", track_name + "_eta", 25, -2.5, 2.5)
    numberOfValidHits = ROOT.TH1F( track_name + "_numberOfValidHits", track_name + "_numberOfValidHits", 25, 0.0, 25.0)

    for event in tree:
        if (eval("event." + branch_variable + "[0]") != 0):
            FillStandardHistos(branch_variable,event,pt,phi,eta,numberOfValidHits)

    print "Histos filled."
    DrawStandardHistos(track_name,pt,phi,eta,numberOfValidHits)
    print "Histos saved.\n"

#-------------------------------------------------------------------------------
# main code

# open root file & tree
filename = "output_gsf_associator_assoc_10000.root"
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
# MakeAssocParaHistos(branch_variable,track_name,assoc_para)
MakeAssocParaHistos("gsf_track","non_assoc_gsf_track",0)
MakeAssocParaHistos("gsf_track","seed_assoc_gsf_track",1)
MakeAssocParaHistos("gsf_track","all_assoc_gsf_track",2)
# MakeQualHistos(branch_variable,track_name)
MakeQualHistos("gsf_track","qual_all_gsf_tracks")
# MakeHistoIfNonZero(branch_variable,track_name)
MakeHistoIfNonZero("seed_assoc_track","seed_assoc_track")
MakeHistoIfNonZero("track_assoc_track","track_assoc_track")

print "All Histos printed."

