
import math
import numpy

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
