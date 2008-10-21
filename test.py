#!/usr/bin/env python

import VibTools

from numpy import *

import pylab as plt

def print_mat (a, fact=1.0, nodes=False) :
    print
    for k in a :
        for l in k :
            print "%6.1f" % (l*fact),

        if nodes :
            n = 0
            for i in range(1,len(k)) :
                if k[i-1]*k[i] < 0.0 :
                    n += 1
            print "nodes: ", n,
        print
    print
 
def test_snf() :
    res = SNFResults()
    res.read()
    res.test_consistency()

def test1_new () :
    res = VibTools.SNFResults()
    res.read()

    # plot IR spectrum

    ir_spectrum = VibTools.VibSpectrum(res.modes, res.get_ir_intensity())

    IRPlot = ir_spectrum.get_plot(1800.0, 1100.0, 0.0, 230.0, spectype='IR')

    IRPlot.scale_range(1100.0, 1400.0, 5.0, boxy=(0.0, 120.0), boxlabel='x 5')
    IRPlot.delete_peaklabels([0,2,9,12])
    
    IRPlot.draw()

    # plot Raman spectrum

    raman_spectrum = VibTools.VibSpectrum(res.modes, res.get_raman_int())

    RamanPlot = raman_spectrum.get_plot(1800.0, 1100.0, 0.0, 18.0, spectype='Raman')

    RamanPlot.delete_peaklabels([0,2,5,11])

    RamanPlot.draw()

    # combine the two

    pl = VibTools.CombinedPlot()
    pl.plot(2,1,[IRPlot, RamanPlot])
    pl.save_plot('spec.pdf')

    plt.show()

def test2() :
    res = VibTools.SNFResults()
    res.read()

    modes = res.modes.get_range(800.0, 1800.0)
    modes = res.modes.get_subset(range(481,500))

    irint = res.get_ir_intensity(modes)
    ri    = res.get_raman_int(modes)
    
    labels = ["%3i %6.1f %10.4f %10.4f    " % (i, modes.freqs[i], irint[i], ri[i]) for i in range(modes.nmodes)]

    modes.print_attype_composition(labels=labels)

    print "IR intensity   : ", irint.sum()
    print "Raman intensity: ", ri.sum()
    
def composition_analysis (modelist, names) :
    res = VibTools.SNFResults()
    res.read()

    for ml, n in zip(modelist, names) :
        print
        print " *** ", n, " *** "
        print

        modes = res.modes.get_subset(ml)
        
        irint = res.get_ir_intensity(modes)
        ri    = res.get_raman_int(modes)
            
        labels = ["%i %6.1f %10.4f %10.4f    " % (i, f, ir, r)
                  for i, f, ir, r in zip(ml, modes.freqs, irint, ri)]

        modes.print_attype_composition(labels)

        print "IR Absorption          : ", irint.sum()
        print "Raman scattering factor: ", ri.sum()

def test3() :
    alpha = True

    if alpha :
        # alpha-helix
        modelist = [[500], range(481,500), [480], range(462,480), [461], range(421,461), range(401,421),
                    [380,381]+range(383,401), [382], [359]+range(361,380), range(340,359)+[360],
                    [319]+range(321,340), [320]]
        names    = ['COOH', 'amide I', 'NH2 bend', 'amide II', 'mixed', 'CH3 asymm', 'CH3 symm',
                    'CA--H (I)', 'COOH', 'CA--H (II)', 'amide III', 'skel stretch', 'COOH']
    else :
        # 310-helix
        modelist = [[500], range(481,500), [480], range(462,480), [461], [421]+range(423,461), [422], range(401,421),
                    [380,381]+range(383,401), [382], range(360,380), range(340,360), [319]+range(321,340), [320]]
        names    = ['COOH', 'amide I', 'NH2 bend', 'amide II', 'mixed', 'CH3 asymm', 'mixed', 'CH3 symm',
                    'CA--H (I)', 'COOH', 'CA--H (II)', 'amide III', 'skel stretch', 'COOH']

    composition_analysis (modelist, names)

def test4() :
    res = VibTools.SNFResults()
    res.read()

    modes = res.modes.get_subset(range(481,500))
    modes.write_g98out()

    modes.print_residue_composition()
    modes.print_attype_composition()


def print_decomposition_analysis(res, modes) :
    ha = VibTools.HugAnalysis(res, 'IR')

    groups, groupnames = res.mol.residue_groups()
    groups2, polyala_keys = res.mol.attype_groups()

    # decompose into residues
    ha.print_group_coupling_matrix(groups, groupnames, modes)
    # decompose into atom types
    ha.print_group_coupling_matrix(groups2, polyala_keys, modes)

    for imode in range(modes.nmodes) :
        print "*** imode %i *** " % (imode+1)
        # decompose into residues
        ha.print_group_coupling_matrix(groups, groupnames, modes, imode)
        # decompose into atom types
        ha.print_group_coupling_matrix(groups2, polyala_keys, modes, imode)

def test5() :
    res = VibTools.SNFResults()
    res.read()

    modes = res.modes.get_subset(range(481,500))

    print_decomposition_analysis(res, modes)

def test6() :
    res = VibTools.SNFResults()
    res.read()

    modes = res.modes.get_subset(range(481,501))
    modes.write_g98out()

    modes.print_residue_composition()
    modes.print_attype_composition()

    print "PM localization "
    
    lv = VibTools.LocVib(modes, 'PM')
    lv.localize(modes)
    lv.sort_by_residue()
    lv.adjust_signs()

    print
    print "coupling matrix: "
    print_mat(lv.get_couplingmat())
    print "transformation matrix to localized modes: "
    print_mat(lv.transmat.transpose(), 100.0, nodes=True)

    lv.locmodes.write_g98out()

    lv.locmodes.print_residue_composition()
    lv.locmodes.print_attype_composition()

def test7() :
    res = VibTools.SNFResults()
    res.read()

    ha = VibTools.HugAnalysis(res, 'ROA', scale=1e4)

    groups, groupnames = res.mol.residue_groups()
    groups2, polyala_keys = res.mol.attype_groups()

    ha.print_group_coupling_matrix(groups, groupnames, res.modes, 491)

    gcm = ha.get_group_coupling_matrix(groups, res.modes, 491)
    fig1 = VibTools.HugAnalysisPlot()
    fig1.plot(gcm, groupnames)
    fig1.draw()

    fig1.save_plot('test1.pdf')

    ha.print_group_coupling_matrix(groups2, polyala_keys, res.modes, 491)

    gcm = ha.get_group_coupling_matrix(groups2, res.modes, 491)
    fig2 = VibTools.HugAnalysisPlot()
    fig2.plot(gcm, polyala_keys)
    fig2.draw()

    fig2.save_plot('test2.pdf')

    plt.show()

def test8() :
    res = VibTools.SNFResults()
    res.read()

    ha = VibTools.HugAnalysis(res, 'ROA', scale=1e4)

    groups, groupnames = res.mol.residue_groups()
    groups2, polyala_keys = res.mol.attype_groups()

    figs = []

    for imode in range(481,500) :
        print " Mode %i " % imode
        ha.print_group_coupling_matrix(groups, groupnames, res.modes, imode)

        gcm = ha.get_group_coupling_matrix(groups, res.modes, imode)
        fig = VibTools.HugAnalysisPlot()
        fig.plot(gcm, groupnames)

        figs.append(fig)

    fig = VibTools.CombinedPlot()
    fig.plot(5,4,figs)
    fig.save_plot('comb.pdf', 'pdf')

    plt.show()

#import profile
#profile.run('test5()')
#profile.run('test5()', 'stats.prof')

#test2()
#test7()

#test1()
test1_new()
