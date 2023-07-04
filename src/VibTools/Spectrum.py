# -*- coding: utf-8 -*-
#
# This file is part of the
# LocVib 1.3 suite of tools for the analysis for vibrational spectra.
# Copyright (C) 2009-2023 by Christoph R. Jacob and others.
#
#    LocVib is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    LocVib is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with LocVib.  If not, see <http://www.gnu.org/licenses/>.
#
# In scientific publications using the LocVib tools, please cite:
#   Ch. R. Jacob, J. Chem. Phys 130 (2009), 084106.
# 
# The most recent version of LocVib is available at
#   http://www.christophjacob.eu/software

import math
import numpy

from . import Plotting
from . import Modes


class VibSpectrum:

    def __init__ (self, freqs, ints) :
        """
        VibSpecctrum constructor.
        """
        if isinstance(freqs, Modes.VibModes) :
            self.freqs = freqs.freqs
        else :
            self.freqs = freqs
        self.ints  = ints

    def get_line_spectrum (self, minfreq, maxfreq, scale=1.0) :
        indices = numpy.where((self.freqs > minfreq) & (self.freqs < maxfreq))
        freqs = self.freqs[indices]
        ints  = self.ints[indices]
        
        x = numpy.zeros((freqs.shape[0]*3+2,))
        x[1:-1:3] = freqs.copy()-0.01
        x[2:-1:3] = freqs.copy()
        x[3:-1:3] = freqs.copy()+0.01

        x[0] = 0.0
        x[-1] = 1.0e6
        
        y = numpy.zeros(x.shape)
        y[2::3] = ints.copy() * scale

        return x, y

    def get_lorentz_spectrum (self, minfreq, maxfreq, halfwidth=15.0) :
        indices = numpy.where((self.freqs > minfreq) & (self.freqs < maxfreq))
        freqs = self.freqs[indices]
        ints  = self.ints[indices]

        x = numpy.arange(minfreq, maxfreq, 0.1)
        y = numpy.zeros(x.shape)
        
        for nu0, intens in zip(freqs, ints) :
            y += ((intens*halfwidth)/(2*math.pi)) / ((x-nu0)**2 + (halfwidth/2)**2)

        return x, y

    def get_gaussian_spectrum (self, minfreq, maxfreq, halfwidth=5.0) :
        indices = numpy.where((self.freqs > minfreq) & (self.freqs < maxfreq))
        freqs = self.freqs[indices]
        ints  = self.ints[indices]

        x = numpy.arange(minfreq, maxfreq, 0.1)
        y = numpy.zeros(x.shape)

        gamma = halfwidth / (2.0*math.sqrt(2.0*math.log(2.0)))
        
        for nu0, intens in zip(freqs, ints) :
            y += (intens/(2.0 * math.pi * gamma**2)) * numpy.exp(-((x-nu0)**2/(2.0*gamma**2)))
            
        return x, y

    def get_rect_spectrum (self, minfreq, maxfreq, bands) :
        x = numpy.arange(minfreq, maxfreq, 0.1)
        y = numpy.zeros(x.shape)

        for b in bands :
            fmin = numpy.min(self.freqs[b])
            fmax = numpy.max(self.freqs[b])
            totint = numpy.sum(self.ints[b])

            indices = numpy.where((x > fmin) & (x < fmax))
            y[indices] = totint
            
        return x, y

    # FIXME: move the following to a VibSpectrumPlot class
    
    def scale_range (self, x, y, minfreq, maxfreq, scale) :
        indices = numpy.where((x > minfreq) & (x < maxfreq))
        y[indices] = y[indices] * scale
        return x, y
        
    def get_band_maxima (self, x, y) :
        maxlist = []
        for n in range(1, len(x)-1) :
            if (y[n] > 0.0) and (y[n] > y[n-1]) and (y[n] > y[n+1]) :
                maxlist.append((x[n], y[n]))
        return maxlist

    def get_band_minima (self, x, y) :
        minlist = self.get_band_maxima(x, -y)
        minlist = [(x,-y) for x,y in minlist]
        return minlist

    def get_plot (self, xmin, xmax, ymin=None, ymax=None, lineshape='Lorentz', hw=None,
                  include_linespec=True, scale_linespec=0.2, spectype='General',
                  label_maxima=True, label_minima=None, boxes=None) :

        if label_minima is None:
            if spectype == 'ROA' :
                label_min=True
            else:
                label_min=False
        else :
            label_min = label_minima

        xmi = min(xmin, xmax)
        xma = max(xmin, xmax)

        spec  = []
        style = []
        lw    = []

        if include_linespec :
            spec.append(self.get_line_spectrum(xmi, xma, scale=scale_linespec))
            style.append('g-')
            lw.append(0.5)
        else :
            spec.append((numpy.array([xmi, xma]), numpy.array([0., 0.])))
            style.append('g-')
            lw.append(0.5)

        if lineshape=='Lorentz' :
            if hw is None :
                spec.append(self.get_lorentz_spectrum(xmi, xma))
            else:
                spec.append(self.get_lorentz_spectrum(xmi, xma, halfwidth=hw))
        elif lineshape=='Gaussian' :
            if hw is None :
                spec.append(self.get_gaussian_spectrum(xmi, xma))
            else:
                spec.append(self.get_gaussian_spectrum(xmi, xma, halfwidth=hw))

        style.append('k-')
        lw.append(2.0)

        xlims = [xmin, xmax]
        ylims = [ymin, ymax]
        if ymin is None:
            ylims[0] = (spec[0][1]).min() * 1.1
        if ymax is None:
            ylims[1] = (spec[0][1]).max() * 1.1

        pl = Plotting.SpectrumPlot()
        pl.plot(spec, style=style, lw=lw, spectype=spectype, xlims=xlims, ylims=ylims)

        if label_maxima :
            maxima = self.get_band_maxima(spec[1][0], spec[1][1]) 
            pl.add_peaklabels(list(zip(maxima, [1]*len(maxima))))
        if label_min :
            minima = self.get_band_minima(spec[1][0], spec[1][1]) 
            pl.add_peaklabels(list(zip(minima, [-1]*len(minima))))

        if boxes is not None :
            bands, names = boxes
            
            for b, n in zip(bands, names) :
                fmin = numpy.min(self.freqs[b])
                fmax = numpy.max(self.freqs[b])
                
                indices = numpy.where((spec[-1][0] > fmin) & (spec[-1][0] < fmax))

                imin = min(0.0, numpy.min(spec[-1][1][indices]))
                imax = max(0.0, numpy.max(spec[-1][1][indices]))

                pl.add_box(fmin, fmax, imin, imax, n)

        return pl

    def get_rect_plot (self, xmin, xmax, ymin=None, ymax=None,
                       bands=None, bandnames=None,
                       include_linespec=True, scale_linespec=0.2,
                       spectype='General') :

        xmi = min(xmin, xmax)
        xma = max(xmin, xmax)

        spec  = []
        style = []
        lw    = []

        if include_linespec :
            spec.append(self.get_line_spectrum(xmi, xma, scale=scale_linespec))
            style.append('g-')
            lw.append(0.5)

        spec.append(self.get_rect_spectrum(xmi, xma, bands))
        style.append('k-')
        lw.append(2.0)

        xlims = [xmin, xmax]
        ylims = [ymin, ymax]
        if ymin is None:
            ylims[0] = (spec[0][1]).min() * 1.1
        if ymax is None:
            ylims[1] = (spec[0][1]).max() * 1.1

        pl = Plotting.SpectrumPlot()
        pl.plot(spec, style=style, lw=lw, spectype=spectype, xlims=xlims, ylims=ylims)

        for b, name in zip(bands, bandnames) :
            fmin = numpy.min(self.freqs[b])
            totint = numpy.sum(self.ints[b])
            pl.labels.append( (fmin, totint, name) )

        return pl

