__author__ = "Andra Stroe"
__version__ = "0.1"

from contextlib import contextmanager
from collections import namedtuple
import gc

import numpy as np
import gleam.matplotlibparams
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import (
    TransformedBbox,
    BboxPatch,
    BboxConnector,
)
from astropy.visualization import quantity_support

import gleam.read_files as rf
import gleam.gaussian_fitting as gf
import gleam.spectra_operations as so

quantity_support()


def mark_inset(parent_axes, inset_axes, loc1a=1, loc1b=1, loc2a=2, loc2b=2, **kwargs):
    """
    Redefined the mark_inset function. Gives freedom to connect any corner of an 
    inset plot to any corner of the parent plot.
    """
    rect = TransformedBbox(inset_axes.viewLim, parent_axes.transData)

    pp = BboxPatch(rect, fill=False, **kwargs)
    parent_axes.add_patch(pp)

    p1 = BboxConnector(inset_axes.bbox, rect, loc1=loc1a, loc2=loc1b, **kwargs)
    inset_axes.add_patch(p1)
    p1.set_clip_on(False)
    p2 = BboxConnector(inset_axes.bbox, rect, loc1=loc2a, loc2=loc2b, **kwargs)
    inset_axes.add_patch(p2)
    p2.set_clip_on(False)

    return pp, p1, p2


def overplot_lines(ax, linelabels, lineloc):
    """
    Overplots emission lines on an already existing plot with axes
    Input:
        ax: matplolib axis
        linelabels: list of latex formatted names for the lines
        lineloc: list of wavelengths of the lines (most likely in Ang)
    Output:
        ax2: returns the new axis in case we want to overplot some more stuff on
             top
    """
    ax2 = ax.twiny()

    for xc in lineloc:
        ax.axvline(x=xc, color="0.8", linestyle="--")

    ax.set_zorder(ax2.get_zorder() + 1)  # put ax in front of ax2
    ax.patch.set_visible(False)  # hide the 'canvas'
    ax2.patch.set_visible(True)
    ax2.axes.get_xaxis().set_ticks(lineloc)
    ax2.axes.xaxis.set_ticklabels(linelabels)
    ax2.set_xlim(ax.get_xlim())
    ax2.xaxis.set_label_position("top")
    ax2.set_xticklabels(ax2.xaxis.get_majorticklabels(), rotation=90)
    return ax2


def overplot_sky(ax, z, sky):
    """
    Overplots sky bands on an already existing plot with axes
    Input:
        ax: matplolib axis
        z: redshift of the source
        sky: sky bands in QTable format
    Output:
        ax: returns the new axis in case we want to overplot some more stuff on
             top
    """
    if sky is None:
        return
    [
        ax.fill_between(
            [
                so.restframe_wl(band["wavelength_min"], z),
                so.restframe_wl(band["wavelength_max"], z),
            ],
            ax.get_ylim()[0],
            ax.get_ylim()[1],
            facecolor="gray",
            alpha=0.3,
        )
        for band in sky
    ]


def plot_spectrum(
    data_path,
    target,
    spectrum,
    spectrum_line,
    spectrum_fit,
    line_names,
    line_latex,
    line_wls,
    fitted_wl,
    inspect,
    cont_width,
    rest_spectral_resolution,
    sky,
):
    """
    Plots spectrum and a fit to a specific line. Overplots the entire line 
    catalogue. 
    Input:
        data_path: folder location of where the plot will be saved
        target: target structure containing info on the target
        spectrum: entire observed spectrum, format from rf.read_spectrum
        spectrum_line: observed spectrum selected around fitted lines
        spectrum_fit: fit to the spectrum (sum of continuum and Gaussians)
        list_names: list of names of the lines that were fit
        line_latex: list of latex formatted line names to overplot
        line_wls: list of the wavelengths of the overplotted lines
        fitted_wl: wavelength(s) of the fitted line
        inspect: do you want to inspect the data dynamically? then it plots the 
                 images with show(); otherwise it prints out the figure into a
                 png file
        plot_fit: do you want to plot only the observed spectrum or overplot the
                  fit to the spectrum
        cont_width: wavelength to the left and right of the line that we will plot
        rest_spectral_resolution: restframed FWHM of the instrument
        sky: sky regions to be masked    
    Output:
        Figure in show() or a saved figure in an external png file
    """
    # Set the title to the plot
    title = (
        f"{target['Sample']}\t".replace("_", "\_")
        + target["Pointing"].replace("_", "\_")
        + f"\t{target['SourceNumber']}\tz={float(target['Redshift']):.3}".replace(
            "_", "\_"
        )
    )

    # Set the basename name of the png outfile
    basename = rf.naming_convention(
        data_path,
        target["Sample"],
        target["SourceNumber"],
        target["Setup"],
        target["Pointing"],
        "linefits",
    )
    for i in range(len(line_names)):
        basename = f"{basename}.{line_names[i].strip()}"

    # Set up figure and axes
    fig = plt.figure(211)
    ax = fig.add_subplot(111)

    # Overplot the observed spectrum
    ax.plot(spectrum["wl_rest"], spectrum["flux_rest"], color="gray")
    # Overplot the observed spectrum around the emission line in a bold colour
    ax.plot(spectrum_line["wl_rest"], spectrum_line["flux_rest"], color="black")

    # plot zoom in on the fitted line
    sub_axes = inset_axes(
        ax,
        width="100%",
        height="100%",
        loc="upper left",
        bbox_to_anchor=(0.5, 1 - 0.35, 0.3, 0.3),
        bbox_transform=ax.transAxes,
    )

    # Select the region around the emission line and overplot the line
    select = (spectrum["wl_rest"] < (np.average(fitted_wl) + cont_width)) & (
        spectrum["wl_rest"] > (np.average(fitted_wl) - cont_width)
    )
    sub_axes.plot(spectrum["wl_rest"][select], spectrum["flux_rest"][select], color="gray")
    ax.set_xlim(
        min(spectrum["wl_rest"]) - cont_width, max(spectrum["wl_rest"]) + cont_width
    )
    ylims = list(ax.get_ylim())
    ylims[0] = max(ylims[0], -2 * np.std(spectrum["flux_rest"].value))
    ax.set_ylim(ylims)
    sub_axes.set_xlim(sub_axes.get_xlim())
    sub_axes.set_ylim(sub_axes.get_ylim())

    # Mark the emission line
    for line in fitted_wl:
        sub_axes.axvline(x=line, color="k", linestyle="-")

    plot_gaussian_fit(
        spectrum["wl_rest"][select], spectrum_fit, ax, rest_spectral_resolution
    )
    plot_gaussian_fit(
        spectrum["wl_rest"][select], spectrum_fit, sub_axes, rest_spectral_resolution
    )

    # Overplot the emission lines of reference
    ax2 = overplot_lines(ax, line_latex, line_wls)

    # Overplot sky bands
    overplot_sky(ax, target["Redshift"], sky)
    overplot_sky(sub_axes, target["Redshift"], sky)

    sub_axes.set_zorder(ax.get_zorder() + 1)  # put ax in front of ax2
    ax.patch.set_visible(False)  # hide the 'canvas'
    sub_axes.patch.set_visible(True)

    # Add title to the zoomed in axis
    sub_axes.set_title(f"{title}")
    ax.set_xlabel(
        r"$\boldsymbol{\lambda}$ "
        + f"({spectrum['wl_rest'].unit.to_string('latex_inline')})",
    )
    ax.set_ylabel(
        r"$F_{\boldsymbol{\lambda}}$"
        + f"({spectrum['flux'].unit.to_string('latex_inline')})",
    )

    # show or save the figure depending on the command line flag
    if inspect:
        fig.show()
        # While not keyboard press, i.e. while it is a mouse press go back into
        # the loop until you actually do press a key
        while not plt.waitforbuttonpress():
            pass
    if not inspect:
        fig.savefig(f"{basename}.png", format="png", bbox_inches="tight")
    fig.clf()
    plt.close(fig)
    gc.collect()


@contextmanager
def overview_plot(
    target, data_path, line_groups, spectrum, cont_width, rest_spectral_resolution, sky
):
    """
    Overview plot of the spectrum of a single target, with zoom-in plots around
    the lines that were fit. The zoom-in plots have the individual Gaussian fits 
    of the lines overplotted. The lines grouped together are also plotted 
    together.
    Input:
        target: properties of the target as read in from rf.read_lof
        data_path: place where the plot will be saved
        line_groups: how many groups of lines are there, i.e. how many zoom-in 
                     plots do we want
        spectrum: spectrum of the source
        cont_width: wavelength to the left and right of the line that we will plot
        rest_spectral_resolution: restframed FWHM of the instrument
        sky: sky regions to be masked
    """

    # Generate the title of the plot from information on the target
    title = (
        f"{target['Sample']}\t".replace("_", "\_")
        + target["Pointing"].replace("_", "\_")
        + f"\t{target['SourceNumber']}\tz={float(target['Redshift']):.3}".replace(
            "_", "\_"
        )
    )
    # Basename that will be used in the savefile name
    basename = rf.naming_convention(
        data_path,
        target["Sample"],
        target["SourceNumber"],
        target["Setup"],
        target["Pointing"],
        "linefits",
    )
    Myxlabel = (
        r"$\boldsymbol{\lambda}$ "
        + f"({spectrum['wl_rest'].unit.to_string('latex_inline')})"
    )
    Myylabel = (
        r"$F_{\boldsymbol{\lambda}}$ "
        + f"({spectrum['flux'].unit.to_string('latex_inline')})"
    )
    # How many zoom-in plots will we have? As many as groups of lines, as the
    # lines grouped together will also be plotted together
    No_plots = len(line_groups)

    # Set figure and add axis
    fig = plt.figure()

    ax = fig.add_subplot(311)

    # Plot the base restframe spectrum to the main axis
    ax.plot(spectrum["wl_rest"], spectrum["flux_rest"], color="k")

    # Set x & y axis limits based on the spectrum
    ax.set_xlim(
        min(spectrum["wl_rest"]) - cont_width, max(spectrum["wl_rest"]) + cont_width
    )
    ylims = list(ax.get_ylim())
    ylims[0] = max(ylims[0], -2 * np.std(spectrum["flux_rest"].value))
    ax.set_ylim(ylims)

    j = 0

    # Function to dynamically add the fits for each emission lines as they are
    # produced/come from a loop
    def plot_line(line, spectrum_fit, rest_spectral_resolution, sky):
        # Counter j to keep track in which subplot we should plot the emission
        # line
        nonlocal j

        # Select a small wavelength range around the emission line where to
        # overplot the fit
        select = (
            spectrum["wl_rest"] < (np.average(line["wavelength"]) + cont_width)
        ) & (spectrum["wl_rest"] > (np.average(line["wavelength"]) - cont_width))

        # Create a zoomed in axis focusing on each line
        # Semi-random x-y ratio that does not match the main plot, but better
        # shows the lines

        Zoom = namedtuple("Zoom", ["x", "y"])
        zoom = Zoom(x=0.9 / No_plots, y=1.5)
        axins = inset_axes(
            ax,
            f"{100*zoom.x}%",
            f"{100*zoom.y}%",
            loc="lower left",
            bbox_to_anchor=(0 + (0 if No_plots < 2 else j / (No_plots - 1)), 1.3, 1, 1),
            bbox_transform=ax.transAxes,
        )

        # Hide axis labels for inset axes
        axins.xaxis.label.set_visible(False)
        if j > 0:
            axins.yaxis.label.set_visible(False)
        # Plot the observed spectrum in the zoomed-in axis
        axins.plot(spectrum["wl_rest"][select], spectrum["flux_rest"][select], color="k")

        # Overplot the gaussian fit to the line in the zoomed-in axis
        plot_gaussian_fit(
            spectrum["wl_rest"][select], spectrum_fit, axins, rest_spectral_resolution
        )
        axins.set_xlim(
            [
                np.average(line["wavelength"]) - cont_width,
                np.average(line["wavelength"]) + cont_width,
            ]
        )
        axins.set_ylim(ylims)

        # Mark the emission lines fitted
        [axins.axvline(x=l["wavelength"], color="gray", linestyle="-") for l in line]
        texts = [
            axins.text(
                l["wavelength"],
                axins.get_ylim()[1] * 1.01,
                f"{l['latex']}",
                ha="center",
                va="bottom",
            )
            for l in line
        ]

        # Adjust for potential label overlap
        if len(texts) > 1:
            adjust_labels(texts, fig, axins, zoom)

        # if it's the first plot, i.e. the left-most plot, label the y axis
        if j == 0:
            axins.set_ylabel(Myylabel)

        # Overplot sky
        overplot_sky(axins, target["Redshift"], sky)

        mark_inset(ax, axins, loc1a=4, loc1b=1, loc2a=3, loc2b=2, fc="none", ec="0.5")

        j = j + 1

    # Return the individual lines plots to the overview plot
    yield plot_line

    # Set title, taking into account inset axes, if they exist
    if j == 0:
        ax.set_title(f"{title}")
    else:
        ax.set_title(f"{title}", y=3)

    # Mark atmosphere absorption
    overplot_sky(ax, target["Redshift"], sky)

    ax.set_ylim(ylims)

    # Set the proper labels
    ax.set_xlabel(Myxlabel)
    ax.set_ylabel(Myylabel)
    # Save the figure under the above-decided name
    plt.savefig(f"{basename}.png", format="png", bbox_inches="tight")
    fig.clf()
    plt.close()
    gc.collect()


def plot_gaussian_fit(wl, spectrum_fit, ax, rest_spectral_resolution):
    """
    Plot the a line fit as a continuum + a Gaussian, whenever the line was 
    detected. Plot a dashed line for upper limits.
    """
    for line_fit in spectrum_fit.lines:
        # Plot detections
        if isinstance(line_fit, gf.Line):
            gauss_part = gf.gauss_function(
                wl,
                line_fit.height.value,
                line_fit.wavelength.value,
                line_fit.sigma.value,
            )
            fit_flux = gauss_part + spectrum_fit.continuum.value
            ax.plot(wl, fit_flux)
        # Plot upper limits
        if isinstance(line_fit, gf.NonDetection) & (
            not isinstance(line_fit, gf.NoCoverage)
        ):
            s = so.fwhm_to_sigma(rest_spectral_resolution)
            gauss_part = gf.gauss_function(
                wl, so.amplitude_to_height(line_fit.amplitude, s), line_fit.restwl, s
            )
            fit_flux = gauss_part + spectrum_fit.continuum.value
            ax.plot(wl, fit_flux, linestyle="--")


def adjust_labels(texts, fig, ax, zoom=None):
    """
    Sometime added text labels for lines overlap because the lines are too close.
    This attempts to move the labels horizontally to allow for better visibility.
    ! Changes are made only for 2 and 3 labels, not more.
    Input:
        texts: lists of text that need to be shuffled around
        fig: parent figure
        ax: parent axis
        zoom: if any zoom was added to the axis with respect to the parent axis.
    Return:
        new positions for text objects.
    """
    transf = ax.transData.inverted()
    unzoom = [1 / z for z in zoom]
    bbox = [
        t.get_window_extent(renderer=fig.canvas.get_renderer())
        .expanded(*(unzoom))
        .transformed(transf)
        for t in texts
    ]

    # Calculate the offset between the right side of the current box and
    # compare it to the the left side of the next box
    difs = [bbox[i].x1 - bbox[i + 1].x0 for i in range(len(bbox) - 1)]

    # Make changes only if there is overlap between text boxes
    if any(dif > 0 for dif in difs):
        if len(texts) == 2:
            offsets = ((-difs[0] / 2), (difs[0] / 2))
        if len(texts) == 3:
            offsets = ((-max(difs[0], 0)), (0), (max(difs[1], 0)))
        for (text, offset) in zip(texts, offsets):
            text.set_x(text.get_position()[0].value + offset)
