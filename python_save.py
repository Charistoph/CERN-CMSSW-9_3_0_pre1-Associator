import os
import sys
import ROOT
from math import sqrt,pi

#-------------------------------------------------------------------------------
# Define Funcitons


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

# call functions


print "Python finished."

