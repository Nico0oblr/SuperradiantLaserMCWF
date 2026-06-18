using PyPlot
using PyCall

function latex_plot_setup(; linewidth_pt=345.0,   # LaTeX \linewidth in pt
    fraction=0.5,           # Fraction of \linewidth for figure width
    font_size_pt=10,         # Font size in pt
    scale=1.0,
    ratio = (sqrt(5) + 1) / 2)               # Scale factor for notebook display
# Convert pt to inches
width_in = linewidth_pt * fraction / 72.27
width_in *= scale
height_in = width_in * ratio # Golden ratio
recommended_linewidth = 0.12 * font_size_pt * scale

# Set matplotlib settings
PyPlot.matplotlib.rc("text.latex", preamble="\\usepackage{physics,bm,amsfonts}")
PyPlot.matplotlib.rc("text", usetex=true)
PyPlot.matplotlib.rc("font", family="serif", size=font_size_pt * scale)
PyPlot.matplotlib.rc("lines", linewidth=recommended_linewidth)
PyPlot.matplotlib.rc("figure", figsize=(width_in, height_in))

PyPlot.matplotlib.rc("xtick.major", size=4)  # length of major ticks
PyPlot.matplotlib.rc("ytick.major", size=4)
PyPlot.matplotlib.rc("xtick.major", pad=4)   # padding between tick and label
PyPlot.matplotlib.rc("ytick.major", pad=4)

axis_linewidth = 0.1 * font_size_pt * scale

PyPlot.matplotlib.rc("axes", linewidth=axis_linewidth)
#PyPlot.matplotlib.rc("xtick", width=axis_linewidth)
#PyPlot.matplotlib.rc("ytick", width=axis_linewidth)
legend_fontsize = 1.0 * font_size_pt * scale
legend_marker_scale = 0.05  # scale down legend line/marker size

PyPlot.matplotlib.rc("legend", fontsize=legend_fontsize)
PyPlot.matplotlib.rc("legend", handlelength=2.0) 



return (width_in, height_in)
end

macro R_str(s)
    s
end