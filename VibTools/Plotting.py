
import math
import numpy
import pylab as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

class HugAnalysisPlot (object) :

    def __init__ (self) :
        self.fig = plt.figure()
        self.ax  = None

    def plot_gcm (self, matrix, groupnames, maxval=None) :
        ngroups = matrix.shape[0]

        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect('equal')

        self.ax.set_frame_on(False)
        rect = Rectangle((0.0,0.0), 1.0, 1.0, ec='w', fc='w', transform=self.fig.transFigure)
        self.ax.add_patch(rect)

        # make the area of the circles proportional to the matrix elements

        if maxval is None :
            maxarea = numpy.abs(matrix).max()
        else:
            maxarea = maxval

        for i in range(ngroups) :
            for j in range(i,ngroups) :

                area = (abs(matrix[i,j])/maxarea) * ((0.4**2)*math.pi)
                rad = math.sqrt(area/math.pi)
                
                if matrix[i,j] < 0.0 :
                    fc = 'w'
                else :
                    fc = 'k'
                
                cir = Circle( (j,i), radius=rad, ec='k', fc=fc)
                self.ax.add_patch(cir)

        line = Line2D((-0.5, ngroups-0.5), (-0.5, -0.5), color='k')
        self.ax.add_line(line)
        line = Line2D((ngroups-0.5, ngroups-0.5), (-0.5, ngroups-0.5), color='k')
        self.ax.add_line(line)

        for i in range(ngroups) :
            # horizontal lines
            line = Line2D((i-0.5, ngroups-0.5), (i+0.5, i+0.5), color='k')
            self.ax.add_line(line)

            # vertical lines
            line = Line2D((i-0.5, i-0.5), (-0.5, i+0.5), color='k')
            self.ax.add_line(line)
                
        self.ax.xaxis.tick_top()
        self.ax.yaxis.tick_right()

        self.ax.set_xticks(range(ngroups))
        self.ax.set_yticks(range(ngroups))

        self.ax.set_xticklabels(groupnames)
        self.ax.set_yticklabels(groupnames)

        for line in self.ax.xaxis.get_ticklines():
            line.set_markersize(0)
        for line in self.ax.yaxis.get_ticklines():
            line.set_markersize(0)

        for label in self.ax.xaxis.get_ticklabels():
            label.set_rotation('vertical')

        self.ax.set_xlim(-0.6, ngroups-0.4)
        self.ax.set_ylim(ngroups-0.4, -0.6)

    def save_plot(self, filename, format='pdf') :
        self.fig.savefig(filename, format=format)

    def __del__ (self) :
        plt.close(self.fig)
