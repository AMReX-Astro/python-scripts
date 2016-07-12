#!/usr/bin/env python

# overplot two different contours from a single 2-d plotfile.  This uses the 
# matplotlib library.
#
# 2012-04-12 M. Zingale

import matplotlib                                                               
matplotlib.use('agg')   

import fsnapshot
import numpy
import pylab
import os
import sys
import getopt
import math
import string
import mpl_toolkits.axes_grid1
import matplotlib.lines

#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile, component1, component2, outFile, 
            minval1, maxval1, minval2, maxval2, 
            ncontours, eps, dpi, 
            xmin, xmax, ymin, ymax,
            label1, label2):

    print             minval1, maxval1, minval2, maxval2

    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    if (outFile == ""):
        outFile = "compare_" + component1 + "_" + component2

        if (not eps):
            outFile += ".png"

        else:
            outFile += ".eps"

    else:
        # make sure the proper extension is used
        if (not eps):
            if (not string.rfind(outFile, ".png") > 0):
                outFile = outFile + ".png"

        else:
            if (not string.rfind(outFile, ".eps") > 0):
                outFile = outFile + ".eps"


    #--------------------------------------------------------------------------
    # read in the data from the plotfile
    #--------------------------------------------------------------------------
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)

    time = fsnapshot.fplotfile_get_time(plotfile)

    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plotfile)

    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*(xmax - xmin)/nx
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*(ymax - ymin)/ny


    if (not nz == -1):
        print "ERROR: 2-d support only"
        sys.exit(2)


    # read in the first components
    data1 = numpy.zeros( (nx, ny), dtype=numpy.float64)
    (data1, err1) = fsnapshot.fplotfile_get_data_2d(plotfile, component1, data1)

    data2 = numpy.zeros( (nx, ny), dtype=numpy.float64)
    (data2, err2) = fsnapshot.fplotfile_get_data_2d(plotfile, component2, data2)

    if not err1 == 0 and err2 == 0:
        sys.exit(2)

    data1 = numpy.transpose(data1)
    data2 = numpy.transpose(data2)

    extent = [xmin, xmax, ymin, ymax]


    #--------------------------------------------------------------------------
    # find data limits, etc.
    #--------------------------------------------------------------------------
    if (not xmin == None):
        extent[0] = xmin

    if (not xmax == None):
        extent[1] = xmax

    if (not ymin == None):
        extent[2] = ymin

    if (not ymax == None):
        extent[3] = ymax


    if (minval1 == None):
        minval1 = numpy.min(data1)

    if (maxval1 == None):
        maxval1 = numpy.max(data1)

    if (minval2 == None):
        minval2 = numpy.min(data2)

    if (maxval2 == None):
        maxval2 = numpy.max(data2)


    levels1 = numpy.linspace(minval1, maxval1, ncontours, endpoint=True)
    levels2 = numpy.linspace(minval2, maxval2, ncontours, endpoint=True)


    print levels1

    #--------------------------------------------------------------------------
    # make the figure
    #--------------------------------------------------------------------------
    cs1 = pylab.contour(x, y, data1, ncontours, colors='b', levels=levels1)
    cs2 = pylab.contour(x, y, data2, ncontours, colors='r', levels=levels2)


    # make the labels -- see http://www.scipy.org/Cookbook/Matplotlib/Legend
    # for this technique
    lines = []
    labels = []

    if (not label1 == None):
        line1 = matplotlib.lines.Line2D(range(10), range(10), linestyle='-', color='b')
        lines.append(line1)
        labels.append(label1)

    if (not label2 == None):
        line2 = matplotlib.lines.Line2D(range(10), range(10), linestyle='-', color='r')
        lines.append(line2)
        labels.append(label2)

        
    pylab.legend(lines, labels)


    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    #pylab.clabel(cs, fontsize=9, inline=1)#, fmt=formatter)

    pylab.axis(extent)

    ax = pylab.gca()
    ax.set_aspect("equal")

    fig1 = ax.get_figure()

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    pylab.xlabel("x")
    pylab.ylabel("y")


    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.5)
    else:
        pylab.savefig(outFile, bbox_inches='tight', pad_inches=0.5)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    ./contour2.py [options] component1 component2 plotfile

    Plot contours of variables "component1" and "component2" from
    the single plotfile.

    Options:

       -o outfile    save the plot to the file outfile

       --m1 value    set the minimum data range for component 1
       --M1 value    set the maximum data range for component 1

       --m2 value    set the minimum data range for component 2
       --M2 value    set the maximum data range for component 2

       -n value      set the number of contours to use

       -x value      x minimum for plot
       -X value      x maximum for plot

       -y value      y minimum for plot
       -Y value      y maximum for plot

       --log         plot the logarithm (base-10) of the data

       --eps         make an EPS plot instead of a PNG

       --dpi value   (PNG only) make the plot with the dpi specified by
                     value

       --label1 str  label for component 1
       --label2 str  label for component 2

    Note: this script requires the fsnapshot.so library, compiled with
    f2py using the GNUmakefile in data_processing/python_plotfile/

    """              
    print usageStr



#==============================================================================
# main
#==============================================================================
if __name__== "__main__":

    outFile = ""
    eps = 0

    minvar1 = None
    maxvar1 = None

    minvar2 = None
    maxvar2 = None

    dpi = 100

    ncontours = 10

    xmin = None
    xmax = None
    ymin = None
    ymax = None

    label1 = None
    label2 = None


    try: opts, next = getopt.getopt(sys.argv[1:], "o:m:n:x:X:y:Y:", 
                                    ["eps","dpi=",
                                     "m1=", "M1=", "m2=", "M2=",
                                     "label1=","label2="])
    except getopt.GetoptError:
        print "invalid calling sequence"
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-o":
            outFile = a

        if o == "--m1":
            try: minvar1 = float(a)
            except ValueError:
                print "invalid value for -m1"
                sys.exit(2)

        if o == "--M1":
            try: maxvar1 = float(a)
            except ValueError:
                print "invalid value for -M1"
                sys.exit(2)

        if o == "--m2":
            try: minvar2 = float(a)
            except ValueError:
                print "invalid value for -m2"
                sys.exit(2)

        if o == "--M2":
            try: maxvar2 = float(a)
            except ValueError:
                print "invalid value for -M2"
                sys.exit(2)

        if o == "-n":
            try: ncontours = int(a)
            except ValueError:
                print "invalid value for -n"
                sys.exit(2)

        if o == "-x":
            try: xmin = float(a)
            except ValueError:
                print "invalid value for -x"
                sys.exit(2)

        if o == "-X":
            try: xmax = float(a)
            except ValueError:
                print "invalid value for -X"
                sys.exit(2)

        if o == "-y":
            try: ymin = float(a)
            except ValueError:
                print "invalid value for -y"
                sys.exit(2)

        if o == "-Y":
            try: ymax = float(a)
            except ValueError:
                print "invalid value for -Y"
                sys.exit(2)

        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                print "invalid value for --dpi"
                sys.exit(2)

        if o == "--label1":
            label1 = a

        if o == "--label2":
            label2 = a



    try: component1 = next[0]
    except IndexError:
        print "ERROR: no component1 specified"
        usage()
        sys.exit(2)    

    try: component2 = next[1]
    except IndexError:
        print "ERROR: no component2 specified"
        usage()
        sys.exit(2)    

    try: plotfile = next[2]
    except IndexError:
        print "ERROR: plotfile not specified"
        usage()
        sys.exit(2)


    do_plot(plotfile, component1, component2, outFile, 
            minvar1, maxvar1, minvar2, maxvar2, 
            ncontours, eps, dpi, 
            xmin, xmax, ymin, ymax,
            label1, label2)
