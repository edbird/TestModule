#!/usr/bin/env python

try:
    import ROOT
    print('imported ROOT from system install')
except ImportError, e:
    import sys
    import os
    cadbrew = None
    for p in os.environ['PATH'].split(':'):
        pp = p.split('/')
        if len(pp[-1]) == 0:
            pp.pop()
        if pp[-2] == 'CadfaelBrew' and pp[-1] == 'bin':
            cadbrew = '/'.join(pp[:-1])
            break
    if cadbrew is None:
        raise ImportError, e
    sys.path.append(os.path.join(cadbrew, 'lib','root'))
    import ROOT
    print('imported ROOT from Cadfaelbrew')

import array
import math
import argparse

c = 299.792458 # speed of light in mm / ns;
m_e = 0.510998910 # electron mass in MeV

#def calc_tof(ene, dEne, tlen):
#    gaminv = m_e / (m_e + ene)
#    beta = math.sqrt(1. - gaminv*gaminv)
#    dBeta = dEne * gaminv * gaminv * gaminv / m_e / beta
#    tof = tlen / beta / c
#    dTof = dEne * tof * gaminv * gaminv * gaminv / beta / beta / m_e
#    return (tof, dTof, beta, dBeta)


def convert(inf):#, outf):
    f = ROOT.TFile(inf)
    t = f.Get("TSD")
    
    # Create storage
    anodic_t0 = ROOT.vector('double')()
    anodic_t1 = ROOT.vector('double')()
    anodic_t2 = ROOT.vector('double')()
    anodic_t3 = ROOT.vector('double')()
    anodic_t4 = ROOT.vector('double')()
    cathodic_t5 = ROOT.vector('double')()
    cathodic_t6 = ROOT.vector('double')()
    
    cell_x = ROOT.vector('double')()
    cell_y = ROOT.vector('double')()
    
    caffe_category = array.array('l',[0])
    
    n_gamma = array.array('i',[0])
    n_positron = array.array('i',[0])
    n_electron = array.array('i',[0])
    n_alpha = array.array('i',[0])
    
    
    # Set storage address
    t.SetBranchAddress('timestamp.anodic_t0', anodic_t0)
    t.SetBranchAddress('timestamp.anodic_t1', anodic_t1)
    t.SetBranchAddress('timestamp.anodic_t2', anodic_t2)
    t.SetBranchAddress('timestamp.anodic_t3', anodic_t3)
    t.SetBranchAddress('timestamp.anodic_t4', anodic_t4)
    t.SetBranchAddress('timestamp.cathodic_t5', cathodic_t5)
    t.SetBranchAddress('timestamp.cathodic_t6', cathodic_t6)
    
    t.SetBranchAddress('timestamp.cell_x', cell_x)
    t.SetBranchAddress('timestamp.cell_y', cell_y)
    
    t.SetBranchAddress('truth.caffe_category', caffe_category)
    
    t.SetBranchAddress('truth.n_gamma', n_gamma)
    t.SetBranchAddress('truth.n_positron', n_positron)
    t.SetBranchAddress('truth.n_electron', n_electron)
    t.SetBranchAddress('truth.n_alpha', n_alpha)
    
    #np = array.array('i',[0])
    #ng = array.array('i',[0])
    #t.SetBranchAddress('particle.nofparticles',np)
    #t.SetBranchAddress('particle.nofgammaonly',ng)
    #tvx = array.array('d',[0.]) # true vertex x
    #tvy = array.array('d',[0.]) # true vertex y
    #tvz = array.array('d',[0.]) # true vertex z
    #t.SetBranchAddress('truth.vertex_x',tvx)
    #t.SetBranchAddress('truth.vertex_y',tvy)
    #t.SetBranchAddress('truth.vertex_z',tvz)
    
    entries = []
    for e in xrange(t.GetEntries()):
        t.GetEntry(e)
        # check number of tracks is 2 and no other gammas (unassociated calorimeter hits)
        #if np[0] != 2 or ng[0] != 0:
        #    continue
        # do this regardless of event topology
        entries.append(e)
    
    #t.SetBranchAddress('particle.charge',Q)

    # Create output histograms
    #c_anodic_t0 = TCanvas('c_anodic_t0', 'c_anodic_t0', 600, 400)
    h_anodic_t0 = ROOT.TH1F('h_anodic_t0', 'h_anodic_t0', 100, 4.76, 4.86) # Re-create timestamp histograms
    # TODO: re-create fit!
    # TODO: change root histogram variable names so my older C++ fit code will work with the new data
     
    
    #c_anodic_t1 = ROOT.TCanvas('c_anodic_t1', 'c_anodic_t1', 600, 400)
    h_anodic_t1 = ROOT.TH1F('h_anodic_t1', 'h_anodic_t1', 100, -10000 + 0.0, 100000 + 50.0)   
        
    for e in entries:
        t.GetEntry(e)
        #print v1[0], v2[0], v1[1], v2[1], ene[0], ene[1]
        #each track starts on foil and ends on calo
        
        for t0 in anodic_t0:
            print(t0)
            h_anodic_t0.Fill(t0)
            
        for t1 in anodic_t1:
            h_anodic_t1.Fill(t1)
            
    
    c_anodic_t0 = ROOT.TCanvas('c_anodic_t0', 'c_anodic_t0', 200, 10, 700, 500)
    h_anodic_t0.Draw();
    c_anodic_t0.SaveAs('c_feast_t0.C')
    c_anodic_t0.SaveAs('c_feast_t0.py')
    c_anodic_t0.SaveAs('c_feast_t0.png')
    c_anodic_t0.SaveAs('c_feast_t0.pdf')

    del c_anodic_t0
    
    c_anodic_t1 = ROOT.TCanvas('c_anodic_t1', 'c_anodic_t1', 200, 10, 700, 500)
    h_anodic_t1.Draw();
    c_anodic_t1.SaveAs('c_feast_t1.C')
    c_anodic_t1.SaveAs('c_feast_t1.py') # TODO
    c_anodic_t1.SaveAs('c_feast_t1.png')
    c_anodic_t1.SaveAs('c_feast_t1.pdf')

    del c_anodic_t1
    
    # TODO: For some reason this code makes the file browser crash when run
    # TODO: Not finished this code, moved to C++ version

    
def main():
    parser = argparse.ArgumentParser(description='conversion')
    parser.add_argument('-i','--input',required=True,help='Input directory name')
    #parser.add_argument('-o','--output',required=True,help='Check if this run can begin')
    args = parser.parse_args()
    convert(args.input)#, args.output)


if __name__ == "__main__":
    main()
