#!/usr/bin/env python

import VibTools

from pylab import *
from numpy import *

def test1 () :
    res = VibTools.SNFResults()
    res.read()

    freqs = res.freqs
    irint = res.get_ir_intensity()
    ri    = res.get_raman_int()
    bi    = res.get_backscattering_int()

    ir_spectrum    = VibTools.VibSpectrum(res.modes, irint)
    raman_spectrum = VibTools.VibSpectrum(res.modes, ri)
    roa_spectrum   = VibTools.VibSpectrum(res.modes, bi)

    # line spectra

    ir_line_x, ir_line_y       = ir_spectrum.get_line_spectrum(1100.0, 1800.0, scale=0.2)
    raman_line_x, raman_line_y = raman_spectrum.get_line_spectrum(1100.0, 1800.0, scale=0.2)
    roa_line_x, roa_line_y     = roa_spectrum.get_line_spectrum(1100.0, 1800.0, scale=0.04)

    # Lorentz spectra

    ir_x, ir_y       = ir_spectrum.get_lorentz_spectrum(1100.0, 1800.0)
    raman_x, raman_y = raman_spectrum.get_lorentz_spectrum(1100.0, 1800.0)
    roa_x, roa_y     = roa_spectrum.get_lorentz_spectrum(1100.0, 1800.0)

    ir_spectrum.scale_range(ir_line_x, ir_line_y, 1100.0, 1400.0, 5.0)
    ir_spectrum.scale_range(ir_x, ir_y, 1100.0, 1400.0, 5.0)
    
    print " *** IR maxima *** "
    for x, y in ir_spectrum.get_band_maxima(ir_x, ir_y) :
        print " %8.2f %8.2f " % (x, y)        

    print " *** Raman maxima *** "
    for x, y in raman_spectrum.get_band_maxima(raman_x, raman_y) :
        print " %8.2f %8.2f " % (x, y)        
        
    figure(1)
    
    subplot(311)
    line1, line2 = plot(ir_line_x, ir_line_y, 'g-', ir_x, ir_y, 'k-')
    line1.set_linewidth(0.5)
    line2.set_linewidth(2.0)

    fill([1400, 1400, 1110, 1110], [100,0,0,100], fill=False)
    text(1400, 105, 'x 5', va='bottom', ha='left')

    axis([1800,1100,0.0,200.0])
    yticks(yticks()[0][1:])
    
    text(1450, 190, 'IR spectrum', va='top', ha='center')
    #xlabel('wavenumber [cm$^{-1}$]')
    ylabel(r'absorption [km/mol]')

    subplot(312)
    line1, line2 = plot(raman_line_x, raman_line_y, 'g-', raman_x, raman_y, 'k-')
    axis([1800,1100,0.0,15.0])
    yticks(yticks()[0][1:])

    line1.set_linewidth(0.5)
    line2.set_linewidth(2.0)

    text(1450, 14, 'Raman spectrum', va='top', ha='center')

    xlabel('wavenumber [cm$^{-1}$]')
    ylabel(r'scattering factor [${\AA}^4$/a.m.u.]')

    subplot(313)
    line1, line2 = plot(roa_line_x, roa_line_y, 'g-', roa_x, roa_y, 'k-')
    axis([1800,1100,-0.015,0.015])
    line1.set_linewidth(0.5)
    line2.set_linewidth(2.0)

    savefig('spec.pdf', type='pdf')
    show()

def test2() :
    res = VibTools.SNFResults()
    res.read()

    modes = res.modes.get_range(800.0, 1800.0)
    modes = res.modes.get_subset(range(481,501))

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
        modelist = [[500], range(481,500), [480], range(462,480), [461], range(421,460), range(401,421),
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

    modes = res.modes.get_subset(range(481,501))
    modes.write_g98out()

    modes.print_residue_composition()
    modes.print_attype_composition()


def print_decomposition_analysis(res, modes) :
    ha = VibTools.HugAnalysis(res, 'Raman')

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

    modes = res.modes.get_subset(range(481,501))

    print_decomposition_analysis(res, modes)


        
#import profile
#profile.run('test5()')
#profile.run('test5()', 'stats.prof')

test2()
test5()

