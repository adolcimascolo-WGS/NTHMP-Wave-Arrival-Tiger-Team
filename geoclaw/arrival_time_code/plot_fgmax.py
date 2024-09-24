"""
Plot fgmax output from GeoClaw run.

"""

import matplotlib.pyplot as plt
import os,sys
import numpy
from clawpack.visclaw import geoplot

#from clawpack.geoclaw import fgmax_tools

sys.path.insert(0, '..')
import fgmax_tools  # local version with changes

outdir='_output'
plotdir='_plots'

if not os.path.isdir(plotdir): 
    os.mkdir(plotdir)

fg = fgmax_tools.FGmaxGrid()
fg.outdir = outdir
data_file = os.path.join(outdir, 'fgmax_grids.data')
fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
fg.read_output()

plt.figure(1,figsize=(8,6))
plt.clf()

clines_zeta = [0.025, 0.3, 0.9] + list(numpy.arange(2,9,2)) #python range (numpy) starts
colors = geoplot.discrete_cmap_1(clines_zeta) # sets colors from geoplot.py

zeta = numpy.where(fg.B>0, fg.h, fg.h+fg.B)   # surface elevation in ocean
plt.contourf(fg.X,fg.Y,zeta,clines_zeta,colors=colors) #contourf (matplotlib - fills in contours with specified levels)
plt.colorbar(label='zeta (m)',shrink=0.7,extend='max')
plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k')  # coastline 0=datum of B   
plt.title("Maximum zeta")

# contour lines of arrival time:
if 1: #0 to exclude arrival time contours
    clines_t = numpy.arange(60,181,5) # This sets the range of countour times, last value is interval 
    arrival_t = fg.arrival_time[0,:,:]/60.  # zeta arrival time in minutes
    clines_t_colors = ([.5,.5,.5],) #.5=gray color 0,0,0=black 1,1,1=white
    con_t = plt.contour(fg.X,fg.Y,arrival_t, clines_t,colors=clines_t_colors) 
    plt.clabel(con_t)
    plt.title("Maximum zeta / eta arrival times")

# fix axes:   x1,x2,y1,y2 = [-122.55, -122.48, 47.605, 47.640]
plt.ticklabel_format(style='plain',useOffset=False)
plt.xticks(rotation=20)
plt.gca().set_aspect(1./numpy.cos(fg.Y.mean()*numpy.pi/180.))
plt.axis([-122.55, -122.48, 47.605, 47.640])

fname = os.path.join(plotdir, "zeta_amplitude_times.png")
plt.savefig(fname)
print("Created ",fname)

# make plots for all arrival times:

atime_labels = ['tfirstETA','tfirstPOS','tfirstNEG','tfirstDRAW',
                'tfirstADVIS','tfirstWARN']
                
for itime in range(6):
    arrival_t = fg.arrival_time[itime,:,:]/60.  # arrival time in minutes
    atime = atime_labels[itime]
    print('range of arrival_t for %s: %g to %g' \
            % (atime,numpy.nanmin(arrival_t),numpy.nanmax(arrival_t)))

    plt.figure(10+itime,figsize=(8,6))
    plt.clf()

    clines_tfirst = numpy.arange(60,181,5) # This sets the range of countour times, last value is interval 
    #colors = numpy.flipud(geoplot.discrete_cmap_1(clines_tfirst)) #flipud reverses color map
    colors = geoplot.discrete_cmap_2(clines_tfirst) #cmap_2 from geoplot.py for arrival times
    plt.contourf(fg.X,fg.Y,arrival_t,clines_tfirst,colors=colors)
    plt.colorbar(shrink=0.7,label='minutes')
    clines_t_colors = ([.5,.5,.5],)
    plt.contour(fg.X,fg.Y,arrival_t,clines_tfirst,
                colors=clines_t_colors, linewidths=0.5)


    #plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k')  # post-seismic coastline
    plt.contour(fg.X,fg.Y,fg.h0,[0.05],colors='g')  # original coastline
    plt.title('arrival time %s' % atime)
    # fix axes:
    plt.ticklabel_format(style='plain',useOffset=False)
    plt.xticks(rotation=20)
    plt.gca().set_aspect(1./numpy.cos(fg.Y.mean()*numpy.pi/180.))
    plt.axis([-122.55, -122.48, 47.605, 47.640])

    fname = os.path.join(plotdir, "arrival_time_%s.png" % atime)
    plt.savefig(fname)
    print("Created ",fname)
