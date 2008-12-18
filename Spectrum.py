
import math
import numpy

import Plotting

class VibSpectrum (object) :

    def __init__ (self, modes, ints) :
        self.freqs = modes.freqs
        self.ints  = ints

    def get_line_spectrum (self, minfreq, maxfreq, scale=1.0) :
        indices = numpy.where((self.freqs > minfreq) & (self.freqs < maxfreq))
        freqs = self.freqs[indices]
        ints  = self.ints[indices]
        
        x = numpy.zeros((freqs.shape[0]*3,))
        x[0::3] = freqs.copy()-0.01
        x[1::3] = freqs.copy()
        x[2::3] = freqs.copy()+0.01

        y = numpy.zeros(x.shape)
        y[1::3] = ints.copy() * scale

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

    # FIXME: move the following to a VibSpectrumPlot class
    
    def scale_range (self, x, y, minfreq, maxfreq, scale) :
        indices = numpy.where((x > minfreq) & (x < maxfreq))
        y[indices] = y[indices] * scale
        return x, y
        
    def get_band_maxima (self, x, y) :
        maxlist = []
        for n in range(1, len(x)-1) :
            if (y[n] > y[n-1]) and (y[n] > y[n+1]) :
                maxlist.append((x[n], y[n]))
        return maxlist

    def get_band_minima (self, x, y) :
        return self.get_band_maxima(x, -y)


    def get_plot (self, xmin, xmax, ymin=None, ymax=None, lineshape='Lorentz', hw=None,
                  include_linespec=True, scale_linespec=0.2, spectype='General',
                  label_maxima=True) :

        xmi = min(xmin, xmax)
        xma = max(xmin, xmax)

        if lineshape=='Lorentz' :
            if hw is None :
                spec = [self.get_lorentz_spectrum(xmi, xma)]
            else:
                spec = [self.get_lorentz_spectrum(xmi, xma, halfwidth=hw)]
        elif lineshape=='Gaussian' :
            if hw is None :
                spec = [self.get_gaussian_spectrum(xmi, xma)]
            else:
                spec = [self.get_gaussian_spectrum(xmi, xma, halfwidth=hw)]

        style = ['k-']
        lw = [2.0]

        if include_linespec :
            spec.append(self.get_line_spectrum(xmi, xma, scale=scale_linespec))
            style.append('g-')
            lw.append(0.5)

        xlims = [xmin, xmax]
        ylims = [ymin, ymax]
        if ymin is None:
            ylims[0] = (spec[0][1]).min() * 1.1
        if ymax is None:
            ylims[1] = (spec[0][1]).max() * 1.1

        pl = Plotting.SpectrumPlot()
        pl.plot(spec, style=style, lw=lw, spectype=spectype, xlims=xlims, ylims=ylims)

        if label_maxima :
            maxima = self.get_band_maxima(spec[0][0], spec[0][1]) 
            pl.add_peaklabels(maxima)

        return pl

