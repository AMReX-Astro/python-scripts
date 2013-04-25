#!/usr/bin/env python

# support doing runtime visualization in 2-d.
#
# a vis.in file is read in that is in the form
#
# [general]
# key = value
#
# [density]
# min = X
# max = Y
# log = [0|1]
# eps = Z
#
# ...
# 
# where each variable to be plotted gets its own block.  We then count
# the number of variables and plot them.
#
# general is a special block meant to control overall plotting
# features.
#
# The "-i file" option can be used to specify a different inputs file
#
# We use the matplotlib Imagegrid to make the plot axes easy to setup.
#
# We look at the aspect ratio of the data to ensure that we use the
# best grid layout.
#
# Note: we always do a plot of 1280x720 pixels (or some multiple thereof).
# as this is 720p HD resolution (good for youtube).
#
# here, if eps > 0, we will clip the data to this value at the lower end.
# This prevents us from taking the log of negative numbers

import matplotlib
matplotlib.use('Agg')   # this is important for batch mode on machines w/o a display
import numpy
import pylab
import os
import sys
import getopt
import ConfigParser
import fsnapshot

import math
import string
from mpl_toolkits.axes_grid1 import ImageGrid



#-----------------------------------------------------------------------------
class variable:

    def __init__(self, name="", minval=None, maxval=None, log=0, eps=-1.0):
        self.name = name
        self.min = minval
        self.max = maxval
        self.log = log
        self.eps = eps
        self.data = None


    def __str__(self):
        if self.min == None:
            minStr = "None"
        else:
            minStr = `self.min`

        if self.max == None:
            maxStr = "None"
        else:
            maxStr = `self.max`

        str = "%s: range = [%s, %s], log = %d" % (self.name, minStr, maxStr, self.log)
        return str


class plotAttr:
    
    def __init__(self, numXlabels=None):
        self.numXlabels = numXlabels


class grid:

    def __init__ (self, xmin=0.0, ymin=0.0, xmax=1.0, ymax=1.0, 
                  dx=0.1, dy=0.1):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.dx = dx
        self.dy = dy


#-----------------------------------------------------------------------------
def parseInfile(inFile):

    vars = []

    pAttr = plotAttr()

    parser=ConfigParser.SafeConfigParser()
    parser.optionxform = str  # case sensitive
    parser.read(inFile)

    if (parser.sections() == []):
        sys.exit("ERROR: no variables defined")

    for section in parser.sections():

        if section == "general":
            # general plot attributes
            for option in parser.options(section):
                
                print "in general: ", option
                if option == "numXlabels":
                    try: value=parser.getint(section,option)
                    except ValueError:
                        sys.exit("invalid numXlabels value")

                    print "setting : ", value
                    pAttr.numXlabels = value
                
            
        else:
            # a variable
            vars.append(variable(section))
        
            for option in parser.options(section):

                if option == "min":
                    try: value=parser.getfloat(section,option)
                    except ValueError:
                        sys.exit("invalid min for %s" % (section))
                        
                    vars[len(vars)-1].min = value

                elif option == "max":
                    try: value=parser.getfloat(section,option)
                    except ValueError:
                        sys.exit("invalid max for %s" % (section))
                        
                    vars[len(vars)-1].max = value

                elif option == "log":
                    try: value=parser.getint(section,option)
                    except ValueError:
                        sys.exit("invalid log for %s" % (section))

                    vars[len(vars)-1].log = value

                elif option == "eps":
                    try: value=parser.getfloat(section,option)
                    except ValueError:
                        sys.exit("invalid eps for %s" % (section))

                    vars[len(vars)-1].eps = value

                else:
                    sys.exit("invalid option for %s" % (section))


        #print vars[len(vars)-1]   # debugging
    return pAttr, vars

    
#-----------------------------------------------------------------------------
def setupAxes(F, aspectRatio, nvar):

    # this is a hack -- the ImageGrid doesn't seem to turn off the 
    # offset text on those axes that don't show the y-axis.  onLeft
    # will hold the axis indices of those axes that have the y-axis
    # on the very left of the figure, and therefore will show the
    # y-axis labels
    onLeft = []

    if (aspectRatio == "h"):

        # for <= 3 variables, do a single column
        # for 4 <= # var <= 6, do two columns

        if (nvar <= 3):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (nvar, 1), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

            # all axes touch the left of the figure
            onLeft = list(range(nvar))

        elif (nvar == 4):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (2, 2), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

            onLeft = [0, 2]

        else:
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (3, 2), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="20%")

            onLeft = [0, 2, 4]

    elif (aspectRatio == "v"):

        # always do 1 row -- just much with the spacings here

        if (nvar <= 4):

            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (1, nvar), direction="row",
                               axes_pad = 0.2 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="3%", cbar_pad="8%")

            onLeft = [0]

        else:

            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (1, nvar), direction="row",
                               axes_pad = 0.2 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

            onLeft = [0]

    else:
        
        # for <= 3 variables, do a single row
        # for 4 <= # var <= 6, do 2 rows. 
        if (nvar <= 3):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (1, nvar), direction="row",
                               axes_pad = 0.2 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="10%")

            onLeft = [0]
            
        elif (nvar == 4):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (2, 2), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

            onLeft = [0, 2]

        else:
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (2, 3), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

            onLeft = [0, 3]

    return axGrid, onLeft


#-----------------------------------------------------------------------------
def doPlot(ax, grd, pAttr, var, yoffset):
    extent = [grd.xmin, grd.xmax, grd.ymin, grd.ymax]

    if var.log:

        if (var.eps > 0):
            # clip the data to prevent logs of negative numbers
            pData = var.data.copy()
            pData[pData < var.eps] = var.eps
            pData = numpy.log10(pData)
        else:
            pData = numpy.log10(var.data)            

        if (not var.min == None): 
            pmin = math.log10(var.min)
        else:
            pmin = None
        if (not var.max == None):
            pmax = math.log10(var.max)
        else:
            pmax = None
    else:
        pData = var.data
        pmin = var.min
        pmax = var.max

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3,3))

    im = ax.imshow(pData, origin="lower", interpolation="nearest",
                   vmin=pmin, vmax=pmax, extent=extent)

    ax.set_title(var.name)

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    print pAttr.numXlabels
    if (not pAttr.numXlabels == None):
        print "here"

        xtickvals = grd.xmin + numpy.arange(pAttr.numXlabels)*(grd.xmax - grd.xmin)/pAttr.numXlabels
        print xtickvals
        ax.set_xticks(xtickvals)


    if (not yoffset):
        ax.yaxis.offsetText.set_visible(False)

    ax.cax.colorbar(im, format=formatter)


#-----------------------------------------------------------------------------
def main(inFile, outFile, double, plotFile):

    # get a list of variable objects that contains the information
    # about what to plot
    pAttr, vars = parseInfile(inFile)

    nvar = len(vars)


    # get and store the grid info
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotFile)
    if (not nz == -1):
        sys.exit("ERROR: cannot read a 3-d dataset")


    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plotFile)

    dx = (xmax - xmin)/nx
    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*dx

    dy = (ymax - ymin)/ny
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*dy

    gridInfo = grid(xmin=xmin, xmax=xmax, 
                    ymin=ymin, ymax=ymax, 
                    dx=dx, dy=dy)


    time = fsnapshot.fplotfile_get_time(plotFile)


    # get the data
    for v in vars:
        data = numpy.zeros( (nx, ny), dtype=numpy.float64 )
        (data, err) = fsnapshot.fplotfile_get_data_2d(plotFile, v.name, data)
        if (not err == 0):
            sys.exit("ERROR: unable to read %s" % (v.name) )

        v.data = numpy.transpose(data)


    # find the aspect ratio:
    #
    # aspectRatio = "h" means horizontal
    #               "v" means vertical
    #               "s" means square (to some degree...)

    if (nx >= 2*ny):
        aspectRatio = "h"
    elif (ny >= 1.5*nx):
        aspectRatio = "v"
    else:
        aspectRatio = "s"



    # setup the figure
    if (double == 1):
        F = pylab.figure(1, (25.6, 14.4)) 
    else:
        F = pylab.figure(1, (12.8, 7.2)) 
    F.clf()

    if (double == 1):
        pylab.rcParams.update({'xtick.labelsize': 20,                              
                               'ytick.labelsize': 20,                              
                               'text.fontsize': 24})                               

        pylab.rc("axes", linewidth=2.0)                                            
        pylab.rc("lines", markeredgewidth=2.0)    
        pylab.rc("font", size=18)
    else:
        pylab.rc("font", size=9)


    # setup the axes
    axGrid, onLeft = setupAxes(F, aspectRatio, nvar)

    # plot the data
    n = 0
    while (n < nvar):
        yoffset = 0
        if n in onLeft:
            yoffset = 1

        doPlot(axGrid[n], gridInfo, pAttr, vars[n], yoffset)
        n += 1


    # 5 variables is a tricky case
    if (nvar == 5 and (aspectRatio == "h" or aspectRatio == "s")):
        # turn off the last axes
        axGrid[5].axis('off')
        axGrid[5].cax.axis('off')


    # write the time   
    print "writing time"
    F.text(0.1, 0.01, "t = %g s" % (time), transform = F.transFigure, color="k")

    # automatically make things look better
    try: F.tight_layout(pad=2.0,w_pad=5.0)  # requires matplotlib >= 1.1
    except:
        pass


    if outFile == None:
        pylab.savefig("%s.png" % (plotFile) )
    else:
        pylab.savefig("%s" % (outFile) )


if __name__ == "__main__":

    # parse the commandline options
    inFile = "vis.in"
    outFile = None
    double = 0

    try: opts, next = getopt.getopt(sys.argv[1:], "i:o:d")
    except getopt.GetoptError:
        sys.exit("ERROR: invalid calling sequence")

    for o, a in opts:
        if o == "-i":
            inFile = a

        if o == "-o":
            outFile = a

        if o == "-d":
            double = 1

    try: plotFile = os.path.normpath(next[0])
    except IndexError:
        sys.exit("ERROR: plotfile not specified")


    main(inFile, outFile, double, plotFile)

