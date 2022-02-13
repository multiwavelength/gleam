import matplotlib

matplotlib.use("TkAgg")
from matplotlib import rc
import matplotlib.pyplot as plt

"""
Some reasonable plotting parameters. Feel free to change as per your taste.
"""

plt.rc("text", usetex=True)
plt.rc(
    "text.latex",
    preamble=r"\usepackage{amsmath} \usepackage{amsmath} \usepackage[mode=text,per-mode=symbol]{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{textgreek} \usepackage{sansmath} \sansmath",
)


###########################################################################
#                                                                         #
# SET PLOTTING PARAMETERS                                                 #
#                                                                         #
# plotting parameters
size_major = 8.0
size_minor = 4.0
thick = 1.0

ms1 = 10
t1 = 0.2
w1 = 0.02
w2 = 0.1
l1 = 0.2


plt.rcParams["axes.linewidth"] = thick
plt.rcParams["xtick.major.size"] = size_major
plt.rcParams["ytick.major.size"] = size_major
plt.rcParams["xtick.minor.size"] = size_minor
plt.rcParams["ytick.minor.size"] = size_minor

# Increase the tick-mark widths as well as the widths of lines
# used to draw marker edges to be consistent with the other figure
# linewidths (defaults are all 0.5)
plt.rcParams["xtick.major.width"] = thick
plt.rcParams["ytick.major.width"] = thick
plt.rcParams["xtick.minor.width"] = thick
plt.rcParams["ytick.minor.width"] = thick
plt.rcParams["lines.markeredgewidth"] = thick
plt.rcParams["lines.linewidth"] = thick
plt.rcParams["lines.antialiased"] = True

plt.rcParams.update({"font.size": 15})
plt.rcParams["axes.titlesize"] = 30

plt.rcParams["figure.figsize"] = 20, 10
plt.rcParams["savefig.dpi"] = 300
#                                                                         #
#                                                                         #
#                                                                         #
###########################################################################
