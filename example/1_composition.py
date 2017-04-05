#!/usr/bin/env python

import VibTools

from numpy import *

def main ():
    res = VibTools.SNFResults(outname='Ala10/snf.out', 
                              restartname='Ala10/restart', 
                              coordfile='Ala10/coord')
    res.read()

    # change here to range of modes which are possibly of interest
    ml = range(140,251)

    modes = res.modes.get_subset(ml)
        
    irint = res.get_ir_intensity(modes)
    ri    = res.get_raman_int(modes)
    bi    = res.get_backscattering_int(modes)

    labels = ["%i %6.1f %10.4f %10.4f %10.4f    " % (i, f, ir, r, b)
              for i, f, ir, r, b in zip(ml, modes.freqs, irint, ri, bi)]

    modes.print_attype2_composition(labels)

main()
