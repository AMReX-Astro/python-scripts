#!/usr/bin/env python

# take two plotfiles and a variable and plot the left half of the first 
# plotfile and the right half of the second plotfile on the same axes.
#
# 2014-03-24 M. Zingale

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
def do_plot(plotfile1, plotfile2, plotfile3, component, outFile, 
            minval, maxval, 
            color1, color2, color3,
            ncontours, eps, dpi, 
            xmin, xmax, ymin, ymax,
            label1, label2, label3):


    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    if (outFile == ""):
        outFile = "compare_" + component

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
    # read in the data from the plotfiles
    #--------------------------------------------------------------------------
    (nx1, ny1, nz1) = fsnapshot.fplotfile_get_size(plotfile1)

    time = fsnapshot.fplotfile_get_time(plotfile1)

    (xmin1, xmax1, ymin1, ymax1, zmin1, zmax1) = \
        fsnapshot.fplotfile_get_limits(plotfile1)

    if (not nz1 == -1):
        print "ERROR: 2-d support only"
        sys.exit(2)


    (nx2, ny2, nz2) = fsnapshot.fplotfile_get_size(plotfile2)

    time = fsnapshot.fplotfile_get_time(plotfile2)

    (xmin2, xmax2, ymin2, ymax2, zmin2, zmax2) = \
        fsnapshot.fplotfile_get_limits(plotfile2)


    if not (nx1 == nx2 and ny1 == ny2 and 
            xmin1 == xmin2 and xmax1 == xmax2 and
            ymin1 == ymin2 and ymax1 == ymax2):
        sys.exit("ERROR: grids don't match")


    if not plotfile3 == None:
        (nx3, ny3, nz3) = fsnapshot.fplotfile_get_size(plotfile3)

        time = fsnapshot.fplotfile_get_time(plotfile3)

        (xmin3, xmax3, ymin3, ymax3, zmin3, zmax3) = \
            fsnapshot.fplotfile_get_limits(plotfile3)

        if not (nx1 == nx3 and ny1 == ny3 and 
                xmin1 == xmin3 and xmax1 == xmax3 and
                ymin1 == ymin3 and ymax1 == ymax3):
            sys.exit("ERROR: grids don't match")
        

    x = xmin1 + numpy.arange( (nx1), dtype=numpy.float64 )*(xmax1 - xmin1)/nx1
    y = ymin1 + numpy.arange( (ny1), dtype=numpy.float64 )*(ymax1 - ymin1)/ny1

    # read in the components
    data1 = numpy.zeros( (nx1, ny1), dtype=numpy.float64)
    (data1, err1) = fsnapshot.fplotfile_get_data_2d(plotfile1, component, data1)

    data2 = numpy.zeros( (nx2, ny2), dtype=numpy.float64)
    (data2, err2) = fsnapshot.fplotfile_get_data_2d(plotfile2, component, data2)

    if not err1 == 0 and err2 == 0:
        sys.exit("ERRORS while reading data")

    data1 = numpy.transpose(data1)
    data2 = numpy.transpose(data2)


    if not plotfile3 == None:
        data3 = numpy.zeros( (nx3, ny3), dtype=numpy.float64)
        (data3, err3) = fsnapshot.fplotfile_get_data_2d(plotfile3, component, data3)

        if not err3 == 0:
            sys.exit("ERRORS while reading data")

        data3 = numpy.transpose(data3)

    extent = [xmin1, xmax1, ymin1, ymax1]


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


    if (minval == None):
        minval = min(numpy.min(data1), numpy.min(data2))

    if (maxval == None):
        maxval = max(numpy.max(data1), numpy.max(data2))

    if not plotfile3 == None:
        if (minval == None):
            minval = min(minval, numpy.min(data3))

        if (maxval == None):
            maxval = max(maxval, numpy.max(data3))

    levels = numpy.linspace(minval, maxval, ncontours, endpoint=True)


    #--------------------------------------------------------------------------
    # make the figure
    #--------------------------------------------------------------------------
    
    # left half of data1
    xx = numpy.array(x[0:nx1/2+1])
    yy = numpy.array(y)
    dd = numpy.array(data1[:,0:nx1/2+1])

    cs1 = pylab.contour(xx, yy, dd, 
                        ncontours, colors=color1, levels=levels)

    xx = numpy.array(x[nx2/2:])
    yy = numpy.array(y)
    dd = numpy.array(data2[:,nx2/2:])

    cs2 = pylab.contour(xx, yy, dd,
                        ncontours, colors=color2, levels=levels)

    if not plotfile3 == None:
        xx = numpy.array(x[0:nx1/2+1])
        yy = numpy.array(y)
        dd = numpy.array(data3[:,0:nx1/2+1])

        cs3 = pylab.contour(xx, yy, dd,
                            ncontours, colors=color3, levels=levels)


    # make the labels -- see http://www.scipy.org/Cookbook/Matplotlib/Legend
    # for this technique
    lines = []
    labels = []

    if (not label1 == None):
        line1 = matplotlib.lines.Line2D(range(10), range(10), linestyle='-', color=color1)
        lines.append(line1)
        labels.append(label1)

    if (not label2 == None):
        line2 = matplotlib.lines.Line2D(range(10), range(10), linestyle='-', color=color2)
        lines.append(line2)
        labels.append(label2)

    if not plotfile3 == None:
        if (not label3 == None):
            line3 = matplotlib.lines.Line2D(range(10), range(10), linestyle='-', color=color3)
            lines.append(line3)
            labels.append(label3)

        
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
    ./contoursplit.py [options] component plotfile1 plotfile2 [plotfile3]

    Plot the left half of plotfile1 and the right half of plotfile2
    on the same contour plot.  If plotfile3 is present, then it is also
    plotted on the left.

    Options:

       -o outfile    save the plot to the file outfile

       -m value      set the minimum data range 
       -M value      set the maximum data range 
    
       -c value      color for left plot
       -C value      color for right plot
       -d            color for the optional plotfile3 on the left

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
       --label3 str  label for component 3

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

    minvar = None
    maxvar = None

    dpi = 100

    ncontours = 10

    xmin = None
    xmax = None
    ymin = None
    ymax = None

    label1 = None
    label2 = None
    label3 = None

    color1 = "b"
    color2 = "r"
    color3 = "k"

    try: opts, next = getopt.getopt(sys.argv[1:], "o:m:M:n:x:X:y:Y:c:C:d:", 
                                    ["eps","dpi=",
                                     "label1=","label2=","label3="])
    except getopt.GetoptError:
        print "invalid calling sequence"
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-o":
            outFile = a

        if o == "-m":
            try: minvar = float(a)
            except ValueError:
                syslexit("invalid value for -m")

        if o == "-M":
            try: maxvar = float(a)
            except ValueError:
                sys.exit("invalid value for -M")

        if o == "-c":
            try: color1 = a
            except ValueError:
                sys.exit("invalid value for -c")

        if o == "-C":
            try: color2 = a
            except ValueError:
                sys.exit("invalid value for -C")

        if o == "-d":
            try: color3 = a
            except ValueError:
                sys.exit("invalid value for -d")

        if o == "-n":
            try: ncontours = int(a)
            except ValueError:
                sys.exit("invalid value for -n")

        if o == "-x":
            try: xmin = float(a)
            except ValueError:
                sys.exit("invalid value for -x")

        if o == "-X":
            try: xmax = float(a)
            except ValueError:
                sys.exit("invalid value for -X")

        if o == "-y":
            try: ymin = float(a)
            except ValueError:
                sys.exit("invalid value for -y")

        if o == "-Y":
            try: ymax = float(a)
            except ValueError:
                sys.exit("invalid value for -Y")

        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                sys.exit("invalid value for --dpi")

        if o == "--label1":
            label1 = a

        if o == "--label2":
            label2 = a

        if o == "--label3":
            label3 = a



    try: component = next[0]
    except IndexError:
        print "ERROR: no component specified"
        usage()
        sys.exit(2)    

    try: plotfile1 = next[1]
    except IndexError:
        print "ERROR: plotfile1 not specified"
        usage()
        sys.exit(2)

    try: plotfile2 = next[2]
    except IndexError:
        print "ERROR: plotfile2 not specified"
        usage()
        sys.exit(2)

    try: plotfile3 = next[3]
    except IndexError:
        plotfile3 = None


    do_plot(plotfile1, plotfile2, plotfile3, component, outFile, 
            minvar, maxvar, 
            color1, color2, color3,
            ncontours, eps, dpi, 
            xmin, xmax, ymin, ymax,
            label1, label2, label3)
