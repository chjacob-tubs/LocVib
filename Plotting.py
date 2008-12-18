
import math
import numpy
import pylab as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.patches import Arrow
from matplotlib.lines import Line2D


class Plot (object) :

    def __init__ (self) :
        self.fig = plt.figure()
        self.ax  = None

    def set_figsize(self, figsize) :
        plt.close(self.fig)
        self.fig = plt.figure(figsize=figsize)

    def plot (self) :
        """ Creates the figure, but does not draw it """
        pass

    def draw (self) :
        """ Creates the figure (for  real). """
        self.ax = self.fig.add_subplot(111)
        self.redraw(self.ax)
        
    def redraw (self, ax) :
        """ Creates the figure inside the given axes. """
        pass

    def save_plot(self, filename, format='pdf') :
        """ Saves the figure to a file. """
        if self.ax is None :
            self.draw()
        self.fig.savefig(filename, format=format)

    def close (self) :
        """ Destroys the figure. """
        plt.close(self.fig)


# CombinedPlot has no redraw methods
# -> it cannot be put into another CombinedPlot

class CombinedPlot (Plot) :

    def __init__ (self) :
        Plot.__init__(self)

    def plot (self, ncols, nrows, plots) :
        plt.close(self.fig)

        w = plots[0].fig.get_figwidth()
        h = plots[0].fig.get_figheight()

        self.fig = plt.figure(figsize=(w*nrows,h*ncols))

        for i, p in enumerate(plots) :
            p.close()
            p.redraw(self.fig.add_subplot(ncols, nrows, i+1))

    def draw (self) :
        pass

    
class HugAnalysisPlot (Plot) :

    def __init__ (self) :
        self.fig = plt.figure(figsize=(8.0,8.0))
        self.ax  = None

    def plot (self, matrix, groupnames, maxval=None, showsum=True) :
        self.matrix = matrix
        self.groupnames = groupnames

        if maxval is None :
            self.maxarea = numpy.abs(matrix).max()
        else:
            self.maxarea = maxval

        self.show_sum = showsum
            
    def redraw (self, ax) :
        ngroups = self.matrix.shape[0]

        ax.set_aspect('equal')

        ax.set_frame_on(False)
        rect = Rectangle((-0.5,-0.5), ngroups+0.5, ngroups+0.5, ec='w', fc='w', fill=True)
        ax.add_patch(rect)

        # make the area of the circles proportional to the matrix elements

        for i in range(ngroups) :
            for j in range(i,ngroups) :

                area = (abs(self.matrix[i,j])/self.maxarea) * ((0.4**2)*math.pi)
                rad = math.sqrt(area/math.pi)
                
                if self.matrix[i,j] < 0.0 :
                    fc = 'w'
                else :
                    fc = 'k'
                
                cir = Circle( (j,i), radius=rad, ec='k', fc=fc)
                ax.add_patch(cir)

        line = Line2D((-0.5, ngroups-0.5), (-0.5, -0.5), color='k')
        ax.add_line(line)
        line = Line2D((ngroups-0.5, ngroups-0.5), (-0.5, ngroups-0.5), color='k')
        ax.add_line(line)

        for i in range(ngroups) :
            # horizontal lines
            line = Line2D((i-0.5, ngroups-0.5), (i+0.5, i+0.5), color='k')
            ax.add_line(line)

            # vertical lines
            line = Line2D((i-0.5, i-0.5), (-0.5, i+0.5), color='k')
            ax.add_line(line)

        if self.show_sum :
            ax.text(0.5, ngroups-1.5, 'total: %5.2f' % self.matrix.sum(), size=28)
            
        ax.xaxis.tick_top()
        ax.yaxis.tick_right()

        ax.set_xticks(range(ngroups))
        ax.set_yticks(range(ngroups))

        ax.set_xticklabels(self.groupnames)
        ax.set_yticklabels(self.groupnames)

        for line in ax.xaxis.get_ticklines():
            line.set_markersize(0)
        for line in ax.yaxis.get_ticklines():
            line.set_markersize(0)

        for label in ax.xaxis.get_ticklabels():
            label.set_rotation('vertical')
            label.set_size(20)

        for label in ax.yaxis.get_ticklabels():
            label.set_size(20)

        ax.set_xlim(-0.6, ngroups-0.4)
        ax.set_ylim(ngroups-0.4, -0.6)


class SpectrumPlot (Plot) :

    def __init__ (self) :
        self.fig = plt.figure(figsize=(8.0,4.0))
        self.ax  = None

        self.spectra = []

        self.has_peaklabels = False
        self.has_scalebox   = False

    def plot (self, spectra, style=None, lw=None, spectype=None, xlims=None, ylims=None) :
        self.spectra = spectra

        if style is None :
            self.spectra_style = ['k-'] * len(self.spectra)
        else:
            self.spectra_style = style

        if lw is None :
            self.spectra_lw = [1.0] * len(self.spectra)
        else:
            self.spectra_lw = lw

        if spectype is None :
            self.type = 'General'
        else:
            self.type = spectype

        if self.type == 'General' :
            self.title  = 'Spectrum'
            self.xlabel = 'wavenumber [cm$^{-1}$]'
            self.ylabel = 'intensity [unknown units]'
        elif self.type == 'IR' :
            self.title  = 'IR Spectrum'
            self.xlabel = 'wavenumber [cm$^{-1}$]'
            self.ylabel = r'absorption [km/mol]'
        elif self.type == 'Raman' :
            self.title  = 'Raman Spectrum'
            self.xlabel = 'wavenumber [cm$^{-1}$]'
            self.ylabel = r'scattering factor [${\AA}^4$/a.m.u.]'
        elif self.type == 'ROA' :
            self.title  = 'ROA Spectrum'
            self.xlabel = 'wavenumber [cm$^{-1}$]'
            self.ylabel = r'intensity difference [${\AA}^4$/a.m.u.]'

        if self.type == 'ROA' :
            self.set_figsize((8.0,6.0))

        if xlims is None :
            self.xlims = (max([ sp[0].max() for sp in spectra]), min([ sp[0].min() for sp in spectra]) ) 
        else:
            self.xlims = xlims

        if ylims is None :
            self.ylims = (min([ sp[1].min() for sp in spectra]), max([ sp[1].max() for sp in spectra]) )
        else:
            self.ylims = ylims

    def scale_range (self, min, max, scale, boxy=None, boxlabel=None) :

        for sp in self.spectra :
            indices = numpy.where((sp[0] > min) & (sp[0] < max))
            sp[1][indices] = sp[1][indices] * scale

        if self.has_peaklabels :
            for p in self.peaks :
                if (p[0] > min) and (p[0] < max) :
                    p[1] = p[1] * scale

        if boxy is not None :
            self.has_scalebox = True
            self.scalebox = (min, max, boxy[0], boxy[1])       # box: (xmin, xmax, ymin, ymax)
            if boxlabel is None :
                self.boxlabel = 'x %f4.1' % scale
            else :
                self.boxlabel = boxlabel

    def add_peaklabels (self, peaks) :
        self.has_peaklabels=True
        self.peaks = []
        for p in peaks :
            self.peaks.append(list(p))

    def delete_peaklabels (self, peaknums) :
        self.peaks = [self.peaks[i] for i in range(len(self.peaks)) if i not in peaknums]

    def redraw (self, ax) :
        ax.set_autoscale_on(False)
        ax.get_figure().subplots_adjust(hspace=0.3)

        self.draw_spectra(ax)
        if self.has_peaklabels :
            self.draw_peaklabels(ax)
        if self.has_scalebox :
            self.draw_scalebox(ax)

    def draw_spectra (self, ax) :
        args = []
        for sp, t in zip(self.spectra, self.spectra_style) :
            args += [sp[0], sp[1], t]

        lines = ax.plot(*args)

        for line, lw in zip(lines, self.spectra_lw) :
            line.set_linewidth(lw)

        ax.set_xlim(self.xlims)
        ax.set_ylim(self.ylims)
        yticks = ax.get_yticks()[1:]
        if yticks[-1] > self.ylims[1] :
            yticks = yticks[:-1]

        ax.set_yticks(yticks)
        
        xtit = 0.5*(self.xlims[0]+self.xlims[1])
        ytit = self.ylims[0] + 0.95*(self.ylims[1]-self.ylims[0])

        ax.text(xtit, ytit, self.title, va='top', ha='center')
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
                
    def draw_peaklabels (self, ax) :
        shift = 0.02 * abs(self.ylims[0]-self.ylims[1])
        for p in self.peaks :
            if p[1] > 0.0 :
                ax.text(p[0],  p[1]+shift, '%4.0f' % p[0], withdash=True,
                        va='center', ha='left', rotation=70,
                        dashlength=10.0, dashdirection=1)
            else:
                raise Exception('Peaklabels for minima not implemented')

    def draw_scalebox (self, ax) :
        x1,x2,y1,y2 = self.scalebox
        ax.fill([x2, x2, x1, x1], [y2, y1, y1, y2], fill=False)

        if (self.xlims[0] < self.xlims[1]) :
            xlab = min(x1,x2)
        else:
            xlab = max(x1,x2)
        ylab = max(y1,y2)
        
        ax.text(xlab, ylab*1.05, self.boxlabel, va='bottom', ha='left')


class ModesPlot (Plot) :

    def __init__ (self) :
        self.fig = plt.figure(figsize=(6.0,10.0))
        self.ax  = None

    def plot (self, modes, freqs, ints, locfreqs, modelabels=None, range=None) :
        self.modes = modes
        self.freqs = freqs
        self.ints = ints
        self.locfreqs = locfreqs

        if range is None :
            self.range = (freqs.max(), freqs.min())
        else:
            self.range = range

        if modelabels is None :
            self.modelabels = range(len(freqs))
        else :
            self.modelabels = modelabels

        self.title  = ''
        self.xlabel = ''
        self.ylabel = 'wavenumber [cm$^{-1}$]'

        self.yticks = None

    def redraw (self, ax) :
        ax.set_frame_on(False)
        ax.set_autoscale_on(False)

        xmin = self.modelabels[0]
        xmax = self.modelabels[-1]

        asize = abs(self.range[1]-self.range[0]) / 20.0
        
        rect = Rectangle((xmin-3.5, self.range[1]), xmax-xmin+7.5, abs(self.range[1]-self.range[0]),
                         ec='w', fc='w', fill=True)
        ax.add_patch(rect)

        line = Line2D((xmin-3.5, xmax+4), (self.range[0], self.range[0]), color='k')
        ax.add_line(line)
        line = Line2D((xmin-3.5, xmin-3.5), (self.range[0], self.range[1]), color='k')
        ax.add_line(line)

        for f in self.locfreqs :
            line = Line2D((xmin-3, xmin-1), (f, f), color='k')
            ax.add_line(line)
        
        for i in range(len(self.freqs)) :
            line = Line2D((xmin-0.5, xmax+0.5), (self.freqs[i], self.freqs[i]), color='k')
            ax.add_line(line)

            ax.text(xmax+1, self.freqs[i], 'int:', va='bottom', ha='left')
            ax.text(xmax+3.6, self.freqs[i], '%7.1f' % self.ints[i], va='bottom', ha='right')
            
            for j in range(len(self.freqs)) :
                arrow = Arrow(self.modelabels[j], self.freqs[i], 0, self.modes[i,j]*asize,
                              width=0.5, fc='k', ec='k')
                ax.add_patch(arrow)
                
        ax.set_xlim((xmin-3.6, xmax+4.1))
        ax.set_ylim((self.range[0]+0.1,self.range[1]-0.1))
        
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        
        ax.set_xticks([xmin-2]+self.modelabels)
        ax.set_xticklabels(['localized\n modes']+self.modelabels)

        if self.yticks is not None:
            ax.set_yticks(self.yticks)
        
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        
        ax.text(0.5*(xmax+xmin), self.range[1]+0.5*asize, self.title,
                va='center', ha='center', size=24)
