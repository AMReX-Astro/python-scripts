#!/usr/bin/env python3

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
# cmap = name
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

from __future__ import print_function

import matplotlib
matplotlib.use('Agg')   # this is important for batch mode on machines w/o a display
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import getopt
import configparser
import fsnapshot

import math
from mpl_toolkits.axes_grid1 import ImageGrid


#-----------------------------------------------------------------------------
class Variable(object):

    def __init__(self, name="", minval=None, maxval=None,
                 log=0, eps=-1.0, cmap=None):

        self.name = name
        self.min = minval
        self.max = maxval
        self.log = log
        self.eps = eps
        self.data = None
        self.cmap = cmap
        self.display_name = None

    def __str__(self):
        if self.min is None:
            min_str = "None"
        else:
            min_str = "{}".format(self.min)

        if self.max is None:
            max_str = "None"
        else:
            max_str = "{}".format(self.max)

        return  "{}: range = [{}, {}], log = {}".format(self.name, min_str, max_str, self.log)


class PlotAttr(object):
    """ attributes for the entire domain, regardless of variable """

    def __init__(self, num_xlabels=None, title=None,
                 xmin=None, xmax=None, ymin=None, ymax=None):
        self.num_xlabels = num_xlabels
        self.title = title

        self.font_size = 9

        # this will be what the user requests, and can be a zoom in on the
        # full domain
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

class Grid(object):
    """ this is a container that holds the information required to interpret
        the data for the plot window.  pa is the PlotAttr object that gives
        the extrema that the user requests.  We compare those to the physical
        extent (which comes in through xmin, xmax, ... """

    def __init__(self, pa,
                 xmin=0.0, ymin=0.0, xmax=1.0, ymax=1.0,
                 nx=-1, ny=-1):

        # these are the actual sizes on disk
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.nx = nx
        self.ny = ny

        # cell-centered grid -- for the whole domain
        self.dx = (xmax - xmin)/nx
        self.x = np.linspace(xmin + 0.5*self.dx, xmax - 0.5*self.dx, nx, endpoint=True)

        self.dy = (ymax - ymin)/ny
        self.y = np.linspace(ymin + 0.5*self.dy, ymax - 0.5*self.dy, ny, endpoint=True)


        # the grid info will refer to the view, not the true size on disk.
        # so we use the values read in from the input file to override
        # the actual extent
        if pa.xmin is not None:
            xmin_use = max(pa.xmin, xmin)
        else:
            xmin_use = xmin

        if pa.xmax is not None:
            xmax_use = min(pa.xmax, xmax)
        else:
            xmax_use = xmax

        if pa.ymin is not None:
            ymin_use = max(pa.ymin, ymin)
        else:
            ymin_use = ymin

        if pa.ymax is not None:
            ymax_use = min(pa.ymax, ymax)
        else:
            ymax_use = ymax

        # these are what we want to plot
        self.pxmin = xmin_use
        self.pxmax = xmax_use
        self.pymin = ymin_use
        self.pymax = ymax_use

        # for the aspect ratio
        self.W = self.pxmax - self.pxmin
        self.H = self.pymax - self.pymin

        # for indexing our uniform data array
        self.ix0 = int((xmin_use - xmin)/self.dx)
        self.iy0 = int((ymin_use - ymin)/self.dy)
        self.ix = int((xmax_use - xmin)/self.dx)
        print("in grid_object: ", self.ix, self.nx)

        self.iy = int((ymax_use - ymin)/self.dy)



#-----------------------------------------------------------------------------
def parse_infile(infile):

    pvars = []

    plt_attr = PlotAttr()

    parser = configparser.SafeConfigParser()
    parser.optionxform = str  # case sensitive
    parser.read(infile)

    if parser.sections() == []:
        sys.exit("ERROR: no variables defined")

    for section in parser.sections():

        if section == "general":
            # general plot attributes
            for option in parser.options(section):

                if option == "num_xlabels":
                    try:
                        value = parser.getint(section, option)
                    except ValueError:
                        sys.exit("invalid num_xlabels value")

                    print("setting : ", value)
                    plt_attr.num_xlabels = value

                if option == "title":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid title value")

                    plt_attr.title = value

                if option == "font_size":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid title value")

                    plt_attr.font_size = int(value)

                if option == "xmin":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid xmin value")

                    plt_attr.xmin = float(value)

                if option == "xmax":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid xmax value")

                    plt_attr.xmax = float(value)

                if option == "ymin":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid ymin value")

                    plt_attr.ymin = float(value)

                if option == "ymax":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid ymax value")

                    plt_attr.ymax = float(value)


        else:
            # a variable
            pvars.append(Variable(section))

            for option in parser.options(section):

                if option == "min":
                    try:
                        value = parser.getfloat(section, option)
                    except ValueError:
                        sys.exit("invalid min for %s" % (section))

                    pvars[-1].min = value

                elif option == "max":
                    try:
                        value = parser.getfloat(section, option)
                    except ValueError:
                        sys.exit("invalid max for %s" % (section))

                    pvars[-1].max = value

                elif option == "log":
                    try:
                        value = parser.getint(section, option)
                    except ValueError:
                        sys.exit("invalid log for %s" % (section))

                    pvars[-1].log = value

                elif option == "eps":
                    try:
                        value = parser.getfloat(section, option)
                    except ValueError:
                        sys.exit("invalid eps for %s" % (section))

                    pvars[-1].eps = value

                elif option == "cmap":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid cmap for %s" % (section))

                    pvars[-1].cmap = value

                elif option == "display_name":
                    try:
                        value = parser.get(section, option)
                    except ValueError:
                        sys.exit("invalid cmap for %s" % (section))

                    pvars[-1].display_name = value

                else:
                    sys.exit("invalid option for %s" % (section))


    return plt_attr, pvars


#-----------------------------------------------------------------------------
def setup_axes(fig, aspect_ratio, nvar):

    # this is a hack -- the ImageGrid doesn't seem to turn off the
    # offset text on those axes that don't show the y-axis.  on_left
    # will hold the axis indices of those axes that have the y-axis
    # on the very left of the figure, and therefore will show the
    # y-axis labels
    on_left = []

    if aspect_ratio == "h":

        # for <= 3 variables, do a single column
        # for 4 <= # var <= 6, do two columns

        if nvar <= 3:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(nvar, 1),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="15%")

            # all axes touch the left of the figure
            on_left = list(range(nvar))

        elif nvar == 4:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(2, 2),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="15%")

            on_left = [0, 2]

        else:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(3, 2),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="20%")

            on_left = [0, 2, 4]

    elif aspect_ratio == "v":

        # always do 1 row -- just much with the spacings here

        if nvar <= 4:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(1, nvar),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="3%",
                            cbar_pad="8%")

            on_left = [0]

        else:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(1, nvar),
                            direction="row",
                            axes_pad=0.2,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="15%")

            on_left = [0]

    else:

        # for <= 3 variables, do a single row
        # for 4 <= # var <= 6, do 2 rows.
        if nvar <= 3:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(1, nvar),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="10%")

            on_left = [0]

        elif nvar == 4:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(2, 2),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="15%")

            on_left = [0, 2]

        else:
            axg = ImageGrid(fig, 111, # similar to subplot(111)
                            nrows_ncols=(2, 3),
                            direction="row",
                            axes_pad=0.5,
                            add_all=True,
                            label_mode="L",
                            share_all=True,
                            cbar_location="bottom",
                            cbar_mode="each",
                            cbar_size="5%",
                            cbar_pad="15%")

            on_left = [0, 3]

    return axg, on_left


#-----------------------------------------------------------------------------
def do_plot(ax, gi, plt_attr, var, yoffset):
    extent = [gi.pxmin, gi.pxmax, gi.pymin, gi.pymax]

    cmap = plt.get_cmap("viridis")

    if var.log:

        if var.eps > 0:
            # clip the data to prevent logs of negative numbers
            pdata = var.data.copy()
            pdata[pdata < var.eps] = var.eps
            pdata = np.log10(pdata)
        else:
            pdata = np.log10(var.data)

        if var.min is not None:
            pmin = math.log10(var.min)
        else:
            pmin = None
        if var.max is not None:
            pmax = math.log10(var.max)
        else:
            pmax = None
    else:
        pdata = var.data
        pmin = var.min
        pmax = var.max

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))

    if var.cmap is not None:
        cmap = var.cmap

    im = ax.imshow(pdata[gi.iy0:gi.iy,gi.ix0:gi.ix],
                   origin="lower", interpolation="nearest",
                   vmin=pmin, vmax=pmax, extent=extent, cmap=plt.get_cmap(cmap))

    if var.display_name is None:
        ax.set_title(var.name)
    else:
        ax.set_title(var.display_name)

    ax.set_xlim(gi.pxmin, gi.pxmax)
    ax.set_ylim(gi.pymin, gi.pymax)

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

    if plt_attr.num_xlabels is not None:
        dx_tick = (gi.pxmax - gi.pxmin)/plt_attr.num_xlabels
        xtickvals = gi.pxmin + np.arange(plt_attr.num_xlabels)*dx_tick
        ax.set_xticks(xtickvals)

    if not yoffset:
        ax.yaxis.offsetText.set_visible(False)

    ax.cax.colorbar(im) #, format=formatter)


#-----------------------------------------------------------------------------
def main(infile, out_file, double, plot_file, eps_out):

    # get a list of variable objects that contains the information
    # about what to plot
    plt_attr, pvars = parse_infile(infile)

    nvar = len(pvars)

    # get and store the grid info
    nx, ny, nz = fsnapshot.fplotfile_get_size(plot_file)
    if not nz == -1:
        sys.exit("ERROR: cannot read a 3-d dataset")

    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plot_file)

    gi = Grid(plt_attr,
              xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
              nx=nx, ny=ny)

    time = fsnapshot.fplotfile_get_time(plot_file)


    # get the data
    for v in pvars:
        data = np.zeros( (nx, ny), dtype=np.float64)
        data, err = fsnapshot.fplotfile_get_data_2d(plot_file, v.name, data)
        if not err == 0:
            sys.exit("ERROR: unable to read {}".format(v.name))

        v.data = np.transpose(data)


    # find the aspect ratio:
    #
    # aspect_ratio = "h" means horizontal
    #                "v" means vertical
    #                "s" means square (to some degree...)

    if gi.W >= 2*gi.H:
        aspect_ratio = "h"
    elif gi.H >= 1.5*gi.W:
        aspect_ratio = "v"
    else:
        aspect_ratio = "s"

    # setup the figure
    if double == 1:
        fig = plt.figure(1, (25.6, 14.4))
    else:
        fig = plt.figure(1, (12.8, 7.2))
    fig.clf()

    if double == 1:
        plt.rcParams.update({'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'text.fontsize': 24})

        plt.rc("axes", linewidth=2.0)
        plt.rc("lines", markeredgewidth=2.0)
        plt.rc("font", size=18)
    else:
        plt.rc("font", size=plt_attr.font_size)


    # setup the axes
    axg, on_left = setup_axes(fig, aspect_ratio, nvar)

    # plot the data
    for n in range(nvar):
        yoffset = 0
        if n in on_left:
            yoffset = 1

        do_plot(axg[n], gi, plt_attr, pvars[n], yoffset)


    # 5 variables is a tricky case (since the grid stores 6)
    if nvar == 5 and (aspect_ratio == "h" or aspect_ratio == "s"):
        # turn off the last axes
        axg[5].axis('off')
        axg[5].cax.axis('off')


    # write the time
    fig.text(0.1, 0.01, "t = %g s" % (time),
             transform=fig.transFigure, color="k")

    # automatically make things look better
    try: fig.tight_layout(pad=2.0, w_pad=5.0)  # requires matplotlib >= 1.1
    except:
        pass

    if plt_attr.title is not None:
        fig.text(0.5, 0.95, plt_attr.title,
                 transform=fig.transFigure, color="k",
                 horizontalalignment="center", fontsize=16)


    if out_file is None:
        if eps_out == 1:
            plt.savefig("%s.eps" % (plot_file))
        else:
            plt.savefig("%s.png" % (plot_file))

    else:
        if eps_out == 1:
            plt.savefig("%s" % (out_file))
        else:
            plt.savefig("%s" % (out_file))



if __name__ == "__main__":

    # parse the commandline options
    infile = "vis.in"
    out_file = None
    double = 0
    eps_out = 0

    try: opts, next_arg = getopt.getopt(sys.argv[1:], "i:o:d", ["eps"])
    except getopt.GetoptError:
        sys.exit("ERROR: invalid calling sequence")

    for o, a in opts:
        if o == "-i":
            infile = a

        if o == "-o":
            out_file = a

        if o == "-d":
            double = 1

        if o == "--eps":
            eps_out = 1

    try: plot_file = os.path.normpath(next_arg[0])
    except IndexError:
        sys.exit("ERROR: plotfile not specified")


    main(infile, out_file, double, plot_file, eps_out)
