
scale_factor = 4/(40/25.4)

small_fig={
     'axes.labelsize': 6,
     'axes.titlesize':0,
     'font.size': 6,
     'legend.fontsize': 6,
     'xtick.labelsize': 6,
     'ytick.labelsize': 6,
     'lines.linewidth' : 1.5/scale_factor,
     'lines.markeredgewidth' : 1/scale_factor,
     'lines.markersize' : 6/scale_factor,
     'font.family' : 'sans-serif',
     'font.sans-serif' : 'Arial',
     'pdf.fonttype' : 42,
     'mathtext.fontset' : 'stix',

     'axes.xmargin' :  .05/scale_factor,  # x margin.  See `axes.Axes.margins`
     'axes.ymargin' :  .05/scale_factor,  # y margin See `axes.Axes.margins`
     'axes.linewidth' :  0.8/scale_factor,     # edge linewidth
     'axes.grid' : True,   # display grid or not
     'axes.titlepad' :  0,#6.0/scale_factor     # pad between axes and title in points
     'axes.labelpad' :  4.0/scale_factor,     # space between label and axis

     ### TICKS
     'xtick.top' : True,   # draw ticks on the top side
     'xtick.bottom' : True,   # draw ticks on the bottom side
     'xtick.major.size' : 3.5/scale_factor,      # major tick size in points
     'xtick.minor.size' : 2/scale_factor,      # minor tick size in points
     'xtick.major.width' : 0.8/scale_factor,    # major tick width in points
     'xtick.minor.width' : 0.6/scale_factor,    # minor tick width in points
     'xtick.major.pad' : 3.5/scale_factor,      # distance to major tick label in points
     'xtick.minor.pad' : 3.4/scale_factor,      # distance to the minor tick label in points
     'xtick.direction' : 'in',    # direction: in, out, or inout
     'xtick.minor.visible' : False,  # visibility of minor ticks on x-axis
     'xtick.major.top' : True,   # draw x axis top major ticks
     'xtick.major.bottom' : True,   # draw x axis bottom major ticks
     'xtick.minor.top' : False,   # draw x axis top minor ticks
     'xtick.minor.bottom' : False,   # draw x axis bottom minor ticks

     'ytick.left' : True,   # draw ticks on the left side
     'ytick.right' : True,  # draw ticks on the right side
     'ytick.major.size' : 3.5/scale_factor,      # major tick size in points
     'ytick.minor.size' : 2/scale_factor,      # minor tick size in points
     'ytick.major.width' : 0.8/scale_factor,    # major tick width in points
     'ytick.minor.width' : 0.6/scale_factor,    # minor tick width in points
     'ytick.minor.pad' : 3.4/scale_factor,      # distance to the minor tick label in points
     'ytick.direction' : 'in',    # direction: in, out, or inout
     'ytick.minor.visible' : False,  # visibility of minor ticks on y-axis
     'ytick.major.left' : True,   # draw y axis left major ticks
     'ytick.major.right' : True,   # draw y axis right major ticks
     'ytick.minor.left' : False,   # draw y axis left minor ticks
     'ytick.minor.right' : False,   # draw y axis right minor ticks

     'figure.facecolor' : 'None',   # figure facecolor; 0.75 is scalar gray
     'figure.edgecolor' : 'None',   # figure edgecolor
     'axes.facecolor' : 'None',
     'grid.linewidth' : 0.8/scale_factor,       # in points
     'grid.alpha' : 0.4,       # transparency, between 0.0 and 1.0
     'backend' : 'pdf'

     }

single_column = {
     'axes.labelsize': 5,
     'axes.titlesize':6,
     'font.size': 5,
     'legend.fontsize': 5,
     'xtick.labelsize': 5,
     'ytick.labelsize': 5,
     'lines.linewidth' : 0.5,
     'lines.markeredgewidth' : 0.5,
     'lines.markersize' : 3,
     'font.family' : 'sans-serif',
     'font.sans-serif' : 'Arial',
     'pdf.fonttype' : 42,
     'mathtext.fontset' : 'stix',

     'axes.linewidth' :  0.5,     # edge linewidth
     'axes.grid' : True,   # display grid or not

     ### TICKS
     'xtick.top' : True,   # draw ticks on the top side
     'xtick.bottom' : True,   # draw ticks on the bottom side
     'xtick.major.size' : 1.4,      # major tick size in points
     'xtick.minor.size' : 0.8,      # minor tick size in points
     'xtick.major.width' : 0.5,    # major tick width in points
     'xtick.minor.width' : 0.3,    # minor tick width in points

     'xtick.direction' : 'in',    # direction: in, out, or inout
     'xtick.minor.visible' : False,  # visibility of minor ticks on x-axis
     'xtick.major.top' : True,   # draw x axis top major ticks
     'xtick.major.bottom' : True,   # draw x axis bottom major ticks
     'xtick.minor.top' : False,   # draw x axis top minor ticks
     'xtick.minor.bottom' : False,   # draw x axis bottom minor ticks

     'ytick.left' : True,   # draw ticks on the left side
     'ytick.right' : True,  # draw ticks on the right side
     'ytick.major.size' : 1.4,      # major tick size in points
     'ytick.minor.size' : 0.8,      # minor tick size in points
     'ytick.major.width' : 0.5,    # major tick width in points
     'ytick.minor.width' : 0.3,    # minor tick width in points
     'ytick.direction' : 'in',    # direction: in, out, or inout
     'ytick.minor.visible' : False,  # visibility of minor ticks on y-axis
     'ytick.major.left' : True,   # draw y axis left major ticks
     'ytick.major.right' : True,   # draw y axis right major ticks
     'ytick.minor.left' : False,   # draw y axis left minor ticks
     'ytick.minor.right' : False,   # draw y axis right minor ticks

     'figure.facecolor' : 'None',   # figure facecolor; 0.75 is scalar gray
     'figure.edgecolor' : 'None',   # figure edgecolor
     'axes.facecolor' : 'None',
     'grid.linewidth' : 0.5,       # in points
     'grid.alpha' : 0.4,       # transparency, between 0.0 and 1.0
     'backend' : 'pdf'
     }