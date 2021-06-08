#####################################################################################
# Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
# certain rights in this software.
#####################################################################################
# NOTICE:
# For five (5) years from 10/21/2019 the United States Government is granted for
# itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and perform publicly and
# display publicly, by or on behalf of the Government. There is provision for the
# possible extension of the term of this license. Subsequent to that period or any
# extension granted, the United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
# data to reproduce, prepare derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so. The specific
# term of the license can be identified by inquiry made to National Technology and
# Engineering Solutions of Sandia, LLC or DOE. NEITHER THE UNITED STATES GOVERNMENT,
# NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING
# SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR
# IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
# USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. Any licensee of this software
# has the obligation and responsibility to abide by the applicable export control laws,
# regulations, and general prohibitions relating to the export of technical data.
# Failure to obtain an export control license or other authority from the Government
# may result in criminal liability under U.S. laws.
# (End of Notice)
####################################################################################


def get_plot_specs():
    """
    Plotting Specifications
    color, linewidth, linestyle, marker, markersize
    """
    plot_specs = {
        "color": "k",
        "linewidth": 0.3,
        # 'linewidth': .8,
        "linestyle": "-",
        "marker": "",
        "markersize": 0,
    }

    #  if self._template._label == 'good':
    #     if self._template._metadata['seismic_phase'] == 'Lg':
    #         clr = 'g'
    #     else:
    #         clr = 'c'
    # elif self._template._label == 'bad':
    #     if self._template._metadata['seismic_phase'] == 'Lg':
    #         clr = 'r'
    #     else:
    #         clr = 'm'
    # else:
    #     clr = 'k'

    return plot_specs


def plot_vertical_line(ax, time, ylims=None, ps=None):
    """

    :param ax:
    :param time: location of vertical line (seconds)
    :param ylims: list with [min_y, max_y]
    :param ps: plot specifications
    :return:
    """
    if ylims is None:
        ylims = ax.get_ylim()

    if ps is None:
        ps = {
            "color": "k",
            "linewidth": 2,
            "linestyle": "--",
            "marker": "",
            "markersize": 0,
        }

    ax.axvline(
        x=time,
        ymin=ylims[0],
        ymax=ylims[1],
        color=ps["color"],
        linewidth=ps["linewidth"],
        linestyle=ps["linestyle"],
    )
