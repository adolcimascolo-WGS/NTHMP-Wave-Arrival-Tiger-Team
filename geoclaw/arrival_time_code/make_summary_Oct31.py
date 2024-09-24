 
"""
Generate table summaries of gauge output from GeoClaw for Guemes Channel

"""

from pylab import *
import os, sys
import copy
import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import dtopotools
import interptools

# SETTING directories where input data is found:

try:
    rootdir = os.environ['CLAW']
except:
    print("*** Need to set environment variable CLAW")
    sys.exit()

print('rootdir = ',rootdir)

outdir = '_output'
plotdir = '_plots'
other_figures_dir    = plotdir + '/_other_figures'
if not os.path.isdir(other_figures_dir):
    os.mkdir(other_figures_dir)
print('Will send plots to other_figures_dir = \n  ', other_figures_dir)

### HERE, will need to change for the MLW run
source_name = 'Cascadia Extended L1 @ MHW'
sealevel = 0.0
###

dtopofile = ('/data1/dtopo/CSZ/CSZ_L1-extended-pmel.tt3')
dtopotype = 3  # format of dtopo file

## PRINTNG directories
print('rootdir is:  ',rootdir)
print('outdir is: ',outdir)
print('plotdir is: ',plotdir)
print('other_figures_dir is: ',other_figures_dir)
print('dtopofile is: ',dtopofile)

fprint_path = os.path.join(other_figures_dir,'gauges_report_timesarrival.txt')
print('Will send output to \n   ',fprint_path)

fprint_file = open(fprint_path,'w')

def fprint(*args):
    # build a string out of all the arguments passed to fprint:
    line = ''
    for arg in args:
        line = line + str(arg)
    fprint_file.write('%s\n' % line)

fprint('\nGAUGE REPORT\n')
fprint('EVENT: Cascadia Extended L1 @ MHW')

figure(400, figsize=(8,8))
figure(500, figsize=(8,8))

gaugeno_list=range(1,36)      #gauges 1 to 35 inclusively  for maximums over entire area

hmax_orig_dry=0.0; hmax_orig_wet=0.0;
hmax_area=0.0; zetamax_area = 0.0; etamax_area=0.0; etamax_area_pquake=0.0;
speedmax_area=0.0; momentummax_area=0.0; mfluxmax_area=0.0;
hmin_area = 9999.0
dhmax_area=0.0

# Read in topography and find the average bottom deformation
dtopo = dtopotools.DTopography()
dtopo.read(dtopofile, dtopotype)


gaugeno_dict={}
for gaugeno in gaugeno_list:
    fprint(' ')
    fprint('Gauge no: ',gaugeno,' was in a gauge array')
    
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    xlong,ylat = gauge.location

    ## Determine the subsidence or uplift at this xlong, ylat
    try:
        dzi_gauge = interptools.interp(xlong,ylat,dtopo.X,dtopo.Y,dtopo.dZ[-1,:,:])
        fprint( '(long,lat) = (',xlong,',',ylat,')' )
        fprint('Subsidence/uplift at the gauge was:  %s m'  % dzi_gauge )
    except:
        dzi_gauge = 0.
        fprint( '(long,lat) = (',xlong,',',ylat,')' )
        fprint('The gauge does not overlap dtopo, so no subsidence/uplift' )
    fprint(' ')

    t = gauge.t / 60.   # convert to minutes
    q = gauge.q
    h = q[0,:]

    if 1: #Set u and v set without changing h
        hcopy = copy.copy(h)
        hcopy = where(hcopy>0.01, hcopy, 1.e6)
        u = q[1,:] / hcopy
        v = q[2,:] / hcopy
    
    s = sqrt(u**2 + v**2)
    momentum = h*s
    mflux = h*s*s
    eta = q[3,:]

    lent = len(t)
    B = eta - h          #B has an entry for every t in the gauge0000x.txt file
                         #Below we find the t when the gauge was turned on.

    ####  The initialization state (at t=0) 
    #Let us say h0, B0, eta0 are the pre-quake values of h, bathymetry, 
    #and eta at the gauge as set initially in GeoClaw. They get set a variety of
    #of ways in GeoClaw, but once set we can retrieve them from the gauge0000x.txt
    #file.

    ####  The post co-seismic state (right after sub and/or uplift happens, around 1.2sec)
    #Let us say h_co, B_co, eta_co are the post-coseismic values of h, bathymetry,
    #and eta at the gauge as computed in GeoClaw.  Note that h_co is just h_0!  That
    #is, h does not change right when sub and/or uplift happens, but the surface elevation eta does.
    #That is, the boat is still sitting on the same amount of water the instant after sub/uplift happens,
    #as the water uplifts or subsides with the boat.
    #B_co = B0 + dzi, and sometimes we refer to B_co as Bpost_quake.

    ###  EASIER THINKING 
    # It is so much easier to think about the amount of water h at a location and envision h
    # changing from the initialization state, or from the co-seismic state, rather than
    # thinking about the gauge being in the water or on the land, eta, zeta, the sealevel of the
    # calculation and all the individual things that can first define h initially and then change it
    # as the computation goes on. This new thinking will then work the same for all gauge locations,
    # no matter whether they were onshore, offshore, in water, dry, etc.  Just think about h.  This
    # will help to define the first positive wave and the first negative wave as well.

    ###  FIRST POSITIVE WAVE  (pos_amt is like 0.05m or some positive value)
    # The first positive wave happens the first time h_0 increases by at least pos_amt.
    # That is, we need to check for the first time when (h - h0) >= pos_amt.
    # (Which components of vector h are bigger than the constant h0 by at least pos_amt?)
    # (Then choose the time of the first component that satisfies this condition.)

    ###  FIRST NEGATIVE WAVE (neg_amt is like 0.05m or some positive value)
    # The first negative wave happens the first time when (h0 - h) >= neg_amt.
    # (The constant h0 is bigger than which components of vector h by at least neg_amt?)
    # (Then choose the time of the first component that satisfies this condition.) 

    ### Lets find B0, h0, eta0.  We need to examine when the gauge was turned on to
    # see what we have to do. If it was turned on at some time greater than 1.2sec, then the
    # first line in the gauge file is for this turn on time, and has the post-quake bathymetry
    # in the first line.  If it was turned on at time t=0, then it has the pre-quake bathymetry
    # in its first line.  In either case, we can use dzi to retrieve the other bathymetry value.
    #
    t_gauge_on = t[0]       #It was turned on at time approx. 1 minute.  t[0] should approx 1 since
                            #seconds was changed to be in minutes at the top of the script. 
                            #t[0] means the time in minutes that is reported in row zero (first row) 
                            #of the gauge file. 
            
    fprint(' This gauge was turned on at following time in minutes: ',t_gauge_on)
    fprint(' ')
    if (t_gauge_on > 0.0):  #turned on at a time after sub/uplift happened
                            #which is greater than 1.2sec but hopefully before effect gets
                            #to the gauge. So if you specified a turn on time say 19min
                            #or so for Guemes, this is what is in the gauge file,
                            #assuming your gauge has not felt the effect yet...
                            #For Geumes, you are fine as there was no sub/uplift
                            #right at your gauge, and hopefully your gauges didn't
                            #change from the initialization state before the gauge
                            #was turned on, due to the subsidence at Anacortes...
     
                            #NEW:  If you turned on the gauge too late--after the
                            #      elevation eta and h started changing, we are hoping
                            #      to catch that below, print a warning message
                            #      and still calculate the correct h0 and eta0
                            #      and max and mins that are in the summary table.

        #Read off the turn on state
        #B_turnon will be B post-quake, but h_turnon and eta_turnon might not be
        #from the state directly after the quake if the gauge was turned on too late
        #but it is this turnon information that we have to work with.
        B_turnon = B[0]
        h_turnon = h[0]
        eta_turnon = eta[0]

        #Find the post-seismic state, B_co, h_co, eta_co, be careful in case gauge was turned on too late
        #
        # 
        #thinking:   What are the possibilities for h_co = h0?  It can be 0, if point started as dry.  This
        #            would happen for most points on land if B0 > sealevel.  It can be
        #            h0=sealevel - B0 since B0 + h0 = sealevel for a lot of the ocean points.  If
        #            it were set by qinit so the initial eta0 is etainit, then etaint = B0 + h0 then 
        #            h0 = etainit - B0.
        #
        #Find the B_co right after subsidence/uplift.  Then use this to find the pre-quake B0.
        B_co = B_turnon         #Correct even if gauge turned on too late
        B0 = B_co - dzi_gauge   #Correct even if gauge turned on too late

        #Find the rest of the post-seismic state h_c0 and eta_co, even if gauge turned on too late

        ######  HERE, if this point's initial eta level, eta0 was set specially to be etainit, then we
        ######        need to set etainit here.  For Guemes, I don't think you had any special etainit
        ######        like for canals, but if this gaugeno was one of those special points, you would set it
        ######        here, in the variable called etainit_level.
        ######
        ######  For this gauge, would put logic here to see if water were to be brought up to sealevel or
        ######  to some other level, called etainit_level.  Let's say etainit_level = sealevel for this gauge
        ######  since we don't have that other special category for Geumes. We could use the gauge's long and
        ######  lat from above or its gaugeno somehow to determine whether etainit_level would be e.g. 2meters
        ######  for a canal, or whether it would be sealevel of the job run.  Also, all our gauges were not
        ######  in regions that were forced to be dry either, which would not permit filling up with water
        ######  to any level.  Since we are inside the gaugeno loop, decisions can be different for each gauge.
        ###
        special_gauge = False  #say we have determined this gauge was not in a canal or had special treatment
        no_force_dry = True    #didn't force any gauges to be dry that would have been filled up to etainit_level.
        if (special_gauge):
            etainit_level = 2.0      #example, never used
        else:
            etainit_level = sealevel
        ######
        ######
        #
        fill_up_condition = (B0 < etainit_level) & no_force_dry 
        if (fill_up_condition):         #fill up to etainit_level initially which is sealevel for Guemes.
            h_co = etainit_level - B0   #recall h_co is the same as h0.
        else:
            h_co = 0.0             #location higher than or equal to sealevel, no water initially
        eta_co = B_co + h_co

        #Calculate the rest of the initialization pre-quake state (backing this out!)
        h0 = h_co                 #h_co doesn't change from h0 right after quake
                                  #from post-coseismic state, and lucky that h_co = h0.

        eta0 = B0+h0              #The initial elevation above MHW of the surface of the water h0 pre-quake

        turnon_diff = abs(eta_turnon - eta_co)
        if (turnon_diff > 0.0001):
            fprint (' WARNING -- GAUGE %5i MAY HAVE BEEN TURNED ON LATE ' %gaugeno)
            fprint (' ')
        

    elif (t_gauge_on == 0.0): 
        #Read off the initialization state, the easy case! No backing out info needed!
        #For nearfield tsunamis, if computational time is not a problem, this is way easier!
        B0 = B[0]
        h0 = h[0]
        eta0 = eta[0]
        B_turnon = B0
        h_turnon = h0
        eta_turnon = eta0

        #Calculate the post co-seismic state, the easy case!
        B_co = B0 + dzi_gauge
        h_co = h[0]
        eta_co = B_co + h_co
        
    B0_long = B0*ones(lent)

    fprint(' Gauge Turn-On State: ')
    fprint(' Bathymetry when the gauge was turned on: ',B_turnon )
    fprint(' Flow depth when the gauge was turned on: ',h_turnon )
    fprint(' Eta when the gauge was turned on: ',eta_turnon)
    fprint(' ')
    fprint(' Initialization State: ')
    fprint(' Initial bathymetry at the gauge before quake was: ',B0 )
    fprint(' Initial depth h0 at the gauge before quake was: ',h0 )
    fprint(' Initial eta  eta0 at the gauge before quake was: ',eta0 )
    fprint(' ')
    fprint(' Co-Seismic State: ')
    fprint(' Post co-seismic bathymetry, B_co, at the gauge right after quake was: ',B_co )
    fprint(' Final bathymetry, B[-1], at the gauge at end of job run was: ',B[-1])
    fprint(' These final bathys should be the same if the gauge is centered on the finest level')
    fprint(' Post co-seismic depth, h_co, at the gauge right after quake was: ',h_co )
    fprint(' Post co-seismic eta, eta_co, at the gauge right after quake was: ',eta_co )
    fprint(' ')

    # We define the shoreline as where the contour of the pre-quake bathymetry
    # B0=0.  We do this even when we are doing a sealevel >0 or a sealevel < 0 job
    # run.  Note that sealevel=0 corresponds to a MHW run and the bathymetry is
    # referenced to MHW so everywhere B0=0 is the MHW shoreline, and that is what
    # we define as the shoreline regardless of the type of jobrun we are doing.
    onshore = (B0 > 0.0)
    offshore = logical_not(onshore)

    # Now, we can calculate zeta at this gauge. 
    # Note that zeta is a vector, one entry for every time in the gauge file.
    # If the pre-quake bathy B0 was offshore, zeta = eta,
    # but if the pre-quake bathy was onshore   zeta = h.
    #
    if (offshore):  #zeta = eta for offshore locations
       zeta = B + h
       zeta0 = B0 + h0
    else:           #zeta = h for onshore locations
       zeta = h
       zeta0= h0
    #

    hmax=h.max()                         
    arg_hmax=h.argmax()
    tmax = t[arg_hmax]
    if (h0 >= hmax):
        hmax = h0
        tmax = 0.0        

    hmin = h.min()     
    arg_hmin=h.argmin()
    tmin = t[arg_hmin]    
    if (h0 <= hmin):
       hmin = h0
       tmin = 0.0

    zetamax=zeta.max()                
    arg_zetamax=zeta.argmax()          
    tzetamax = t[arg_zetamax]          
    if (zeta0 >= zetamax):
        zetamax = zeta0
        tzetamax = 0.0

    zetamin=zeta.min()
    arg_zetamin=zeta.argmin()
    tzetamin = t[arg_zetamin]
    if (zeta0 <= zetamin):
        zetamin = zeta0
        tzetamin = 0.0

    dhmax=(zetamax-dzi_gauge)                         
    arg_dhmax=dhmax.argmax()
    if (h0 >= hmax):
        dhmax = h0

    ###  FIRST POSITIVE WAVE
    pos_amt = .1016                               #4 centimeter increase from h0 flags arrival
    hdeep_index = where( (h-h0) > 0.025)[0]      #0.025 is ~1 inches. Works when h0 is positive as well as 0.
                                                #which we need for a sealevel > 0 runs or canals initialized
                                                #with h0>0, for example.
                                                #h is the vector in the gauge file, so starts at turnon time,
                                                #h0 is correct, but we may have exceeded 0.05 above h0 before
                                                #turnon time, and will never know, but will get a hint if this
                                                #gets flagged at turnon time.  Can also look at the values of
                                                #the states that are printed out above the summary table for
                                                #each gauge.
    if len(hdeep_index) >= 1:
        arr_index = hdeep_index[0]
        tfirstpos = t[arr_index]
        tfirstpos = str(round(tfirstpos,1)).rjust(5)  #converted the float to a string with one decimal digit
    else:
        tfirstpos = "  n/a"                       

    ### FIRST NEGATIVE WAVE
    neg_amt = .1016                               #4 centimeter decrease from h0 flags arrival
    hdeep_index = where( (h0-h) > 0.025)[0]      #0.025 is ~1 inches. Might want this to be more. 
                                                #Again, h is the vector in the gauge file that starts at turnon.
                                                #h0 is correct, but h0 really have exceeded h before turnon time.
                                                #Good to check the printout of the initial states for each gauge.
    if len(hdeep_index) >= 1:
        arr_index = hdeep_index[0]
        tfirstneg = t[arr_index]
        tfirstneg = str(round(tfirstneg,1)).rjust(5)     #converted the float to a string with one decimal digit
    else:
        tfirstneg = "  n/a"

    ### FIRST SIGNIFICANT DRAWDOWN
    draw_amt = 0.3048                               #1 foot?
    hdeep_index = where( (h0-h) > 0.3048)[0]     
                                                
                                                
                                               
    if len(hdeep_index) >= 1:
        arr_index = hdeep_index[0]
        tfirstdraw = t[arr_index]
        tfirstdraw = str(round(tfirstdraw,1)).rjust(5)     #converted the float to a string with one decimal digit
    else:
        tfirstdraw = "  n/a"    


    ### FIRST ADVISORY-LEVEL WAVE
    advis_amt = 0.3048                               #1 foot represents advisory threshold
    hdeep_index = where( (h-h0) > 0.3048)[0]     
                                                
                                                
                                               
    if len(hdeep_index) >= abs(1):
        arr_index = hdeep_index[0]
        tfirstadvis = t[arr_index]
        tfirstadvis = str(round(tfirstadvis,1)).rjust(5)     #converted the float to a string with one decimal digit
    else:
        tfirstadvis = "  n/a"  

    ### FIRST WARNING-LEVEL WAVE
    warn_amt = 0.9144                               #3 feet represents warning threshold
    hdeep_index = where( (h-h0) > 0.9144)[0]     
                                                
                                                
                                               
    if len(hdeep_index) >= abs(1):
        arr_index = hdeep_index[0]
        tfirstwarn = t[arr_index]
        tfirstwarn = str(round(tfirstwarn,1)).rjust(5)     #converted the float to a string with one decimal digit
    else:
        tfirstwarn = "  n/a"  
                     

    etamax=eta.max()
    if (eta0 > etamax):
        etamax = eta0
    etamax_pquake = hmax + B_co                 #could have used hmax + B_co, hopefully the same for centered gauges on
                                                #the finest level.  Let's change to use B_c0 instead of B[-1].
    speedmax=s.max()
    momentummax=momentum.max()
    mfluxmax=mflux.max()
    fprint(' Tsunami Data: ')
    fprint(' Maximum value of zeta at the gauge was: ',zetamax )
    fprint(' Maximum value of h, the flow depth, at the gauge was: ',hmax )
    fprint(' Maximum value of dh, the change in water depth, at the gauge was: ',dhmax )
    fprint(' Minimum value of zeta at the gauge was: ',zetamin )
    fprint(' Minimum value of h, the water depth, at the gauge was: ',hmin )
    fprint(' Maximum value of eta, height above MHW, at the gauge was: ',etamax )
    fprint(' Maximum value of eta postquake (hmax + B), height above MHW, at the gauge was: ',\
             etamax_pquake)
    fprint(' Maximum speed at the gauge was: ',speedmax )
    fprint(' Maximum momentum at the gauge was: ',momentummax )
    fprint(' Maximum momentum flux at the gauge was: ',mfluxmax )
    fprint(' ')
    fprint(' Timing Thresholds (in meters): ')
    fprint(' pos_amt and neg_amt were: ',pos_amt,', ',neg_amt)
    fprint(' draw_amt was: ',draw_amt)
    fprint(' advis_amt was: ',advis_amt)
    fprint(' warn_amt was: ',warn_amt)
    fprint(' ')
    fprint(' Time of first positive arrival ((h-h0)>pos_amt) tfirstPOS (in minutes): ',tfirstpos )
    fprint(' Time of first negative arrival ((h0-h)>neg_amt) tfirstNEG (in minutes): ',tfirstneg )
    fprint(' Time of first significant drawdown ((h0-h)>draw_amt) tfirstDRAW (in minutes): ',tfirstdraw )
    fprint(' Time of first advisory-level arrival ((h-h0)>advis_amt) tfirstADVIS (in minutes): ',tfirstadvis )
    fprint(' Time of first warning-level arrival ((h-h0)>warn_amt) tfirstWARN (in minutes): ',tfirstwarn ) 
    fprint(' ')
    fprint(' Maximum value of zeta occurred at time (in minutes): ',tzetamax )
    fprint(' Maximum value of h occurred at time (in minutes): ',tmax )
    fprint(' Minimum value of zeta occurred at time (in minutes): ',tzetamin )
    fprint(' Minimum value of h occurred at time (in minutes): ',tmin )


    if 0:
        figure(400)
        clf()
        subplot(311)
        plot(t, h, 'b')
        xlabel('')
        ylabel('Flow depth (m)')
        title('Gauge %i long=%s,lat=%s,dzi=%5.2f,max h=%5.2f,max eta=%5.2f' %(gaugeno,\
                    xlong,ylat,dzi_gauge,hmax,etamax))

        subplot(312)
        plot(t, zeta, 'b')
        xlabel('')
        ylabel('Zeta (m): h for B0>0; eta for B0<0')
        title('')

        ## eta is always the height of water above the fixed datum called MHW.
        subplot(313)
        plot(t, eta, 'b')
        plot(t, B0_long,'y')
        plot(t, B,'g')
        plot([t[0],t[0]],[B0,B[-1]],'g--')
        xlabel('')
        ylabel('Eta (m): h+B, B0(Y), B(G)')
        title('')

        tight_layout()

        fdir  = other_figures_dir
        fname_gauge = 'Gauge%s.png' % str(gaugeno).zfill(5)
        fname = source_name + '_' + 'part1_' + fname_gauge
        fwhere = fdir + '/' + fname
        savefig(fwhere, bbox_inches='tight')
        fprint('Created %s' % fname )

    if 0: #Choose 0 if already have the plots
        figure(500)
        clf()
        subplot(311)
        plot(t, s, 'b')
        xlabel('')
        ylabel('speed (m/s)')
        title('Gauge %i long=%s,lat=%s,s_max=%5.2f,hs_max=%5.2f,hss_max=%5.2f' %(gaugeno,\
                    xlong,ylat,speedmax,momentummax,mfluxmax))

        subplot(312)
        plot(t, momentum, 'b')
        ylabel('momentum (m^2 / s)')
        title('')

        subplot(313)
        plot(t, mflux, 'b')
        xlabel('time (Minutes after earthquake)')
        ylabel('momentum flux (m^3 / s^2)')
        title('')

        tight_layout()

        fdir  = other_figures_dir
        fname_gauge = 'Gauge%s.png' % str(gaugeno).zfill(5)
        fname = source_name + '_'+ 'part2_' + fname_gauge
        fwhere = fdir + '/'+ fname
        savefig(fwhere, bbox_inches='tight')
        fprint('Created %s' % fname )

    ## Check to see if the maximums over this area were
    ## changed by this gauge.
    if ((hmax > hmax_orig_dry) & (B0 > 0)):
        hmax_orig_dry = hmax
        hmax_orig_dry_gauge = gaugeno      #update the gaugeno with the maximum
    if ((hmax > hmax_orig_wet) & (B0 <= 0)):
        hmax_orig_wet = hmax
        hmax_orig_wet_gauge = gaugeno      #update the gaugeno with the maximum
    if (hmax > hmax_area):
        hmax_area = hmax
        hmax_area_gauge = gaugeno          #update the gaugeno with the maximum
    if (dhmax > dhmax_area):
        dhmax_area = dhmax
        dhmax_area_gauge = gaugeno          #update the gaugeno with the maximum
    if (zetamax > zetamax_area):
        zetamax_area = zetamax
        zetamax_area_gauge = gaugeno       #update the gaugeno with the maximum
    if (etamax > etamax_area):
        etamax_area = etamax
        etamax_area_gauge = gaugeno             #update the gaugeno with the maximum
    if (etamax_pquake > etamax_area_pquake):
        etamax_area_pquake = etamax_pquake
        etamax_area_pquake_gauge = gaugeno #update the gaugeno with the maximum
    if (speedmax > speedmax_area):
        speedmax_area = speedmax
        speedmax_area_gauge = gaugeno      #update the gaugeno with the maximum
    if (momentummax > momentummax_area):
        momentummax_area = momentummax
        momentummax_area_gauge = gaugeno   #update the gaugeno with the maximum
    if (mfluxmax > mfluxmax_area):
        mfluxmax_area = mfluxmax
        mfluxmax_area_gauge = gaugeno           #update the gaugeno with the maximum

    #Find out if this particular hmin is smaller than hmin_area. Only
    #do it for offshore_gauges. (We did calculate hmin for onshore
    #gauges as well, but the min of those is probably 0.)
    if ((hmin < hmin_area) & offshore):
        hmin_area = hmin
        hmin_area_gauge = gaugeno          #update the gaugeno with hmin the smallest
   

 ## Save the info for this gauge in the dictionary below for later printing
    value_dict={'B0': B0, 'B': B_co, 'max_h': hmax, 'max_dh':dhmax, 'min_h': hmin,'max_zeta': zetamax, 'max_eta': etamax,\
            'max_eta_pquake': etamax_pquake, 'max_speed': speedmax, 'max_momentum': momentummax,\
            'max_mflux': mfluxmax, 'dzi': dzi_gauge, 'tmax': tmax, 'tmin': tmin, 'tfirstpos': tfirstpos, 'tfirstneg': tfirstneg,\
            'tfirstdraw': tfirstdraw, 'tfirstadvis': tfirstadvis, 'tfirstwarn': tfirstwarn}
    gaugeno_dict[gaugeno]=value_dict

    if 0:
            figure(400)
            clf()
            subplot(311)
            plot(t, h, 'b')
            xlabel('')
            ylabel('Flow depth (m)')
            title('Gauge %i long=%s,lat=%s,dzi=%5.2f,max h=%5.2f,max eta=%5.2f' %(gaugeno,\
                        xlong,ylat,dzi_gauge,hmax,etamax))

###HERE, the loop over gauges has now completed, so we know the max over all of them and the
####     minimum of hmin over the offshore ones, so print these variables we saved as we went
####     through the loop over gauges

### HERE, I am now printing out the gaugeno that gave the maximum (or the minimum for hmin).
###       You will see the variable that I saved as the gaugeno loop was traversed.  If the current
###       gaugeno created a bigger value than all others before it, that gaugeno was recorded, so at
###       the end of the loop the variable has the gaugeno that is being printed out below.
### Print the maximums over the area encompassed by all the gauges used
fprint(' ')
fprint('   MAXIMUM VALUES OVER ALL 35 GAUGES ' )
fprint(' ' )
fprint(' Maximum value of h (flow depth) over all gauges: %8.2f occurred at gauge: %3i '\
       %(hmax_area,hmax_area_gauge) ) 
#fprint(' Maximum value of h (flow depth) over all onshore gauges (B0>0): %8.2f occurred at gauge: %3i '\
 #      %(hmax_orig_dry,hmax_orig_dry_gauge) ) 
fprint(' Maximum value of h (flow depth) over all offshore gauges (B0<=0): %8.2f occurred at gauge: %3i '\
       %(hmax_orig_wet,hmax_orig_wet_gauge) )
fprint(' Maximum value of dh (change in flow depth) over all gauges: %8.2f occurred at gauge: %3i '\
       %(dhmax_area,dhmax_area_gauge) )  
fprint(' Maximum value of zeta, same as flow depth if onshore: %8.2f occurred at gauge: %3i '\
       %(zetamax_area,zetamax_area_gauge) ) 
fprint(' Maximum value of eta, height above MHW: %8.2f occurred at gauge: %3i '\
       %(etamax_area,etamax_area_gauge) ) 
fprint(' Maximum value of eta post quake, height above MHW: %8.2f occurred at gauge: %3i '\
       %(etamax_area_pquake,etamax_area_pquake_gauge) ) 

fprint(' Minimum value of h (water depth) over offshore gauges: %8.2f occurred at gauge: %3i '\
         %(hmin_area,hmin_area_gauge) )
#####

fprint(' Maximum speed: %8.2f  occurred at gauge: %3i ' %(speedmax_area,speedmax_area_gauge) )
fprint(' Maximum momentum: %8.2f  occurred at gauge: %3i ' %(momentummax_area,momentummax_area_gauge) )
fprint(' Maximum momentum flux: %8.2f  occurred at gauge: %3i ' %(mfluxmax_area,mfluxmax_area_gauge) )
fprint(' ')

### Print summary table all gauges
fprint('               SUMMARY FOR EACH OF THESE GAUGES                ' )
fprint('Gauge    B0       B         dzi     max       min     max       max      max      max     max      max     tmax     tmin  tfirstPOS tfirstNEG tfirstDRAW tfirstADVIS tfirstWARN' )
fprint('                                     h         h      zeta      dh       eta       s      hs       hss                                       ' )
fprint('                                                                      post-quake                                                             ' )
for key in gaugeno_dict: 
    value_dict = gaugeno_dict[key]

    #fprint('%5i %8.3f %8.3f %8.3f %5.2f %5.2f %6.2f %6.2f  %6.2f  %6.2f  %6.2f  %8.2f %8.1f %8.1f %8.1f %8.1f %8.1f' %(key,value_dict['B0'],\

    fprint('%3i %8.2f %8.2f  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f   %s   %s      %s      %s      %s' %(key,value_dict['B0'],\
           value_dict['B'],value_dict['dzi'],value_dict['max_h'],value_dict['min_h'],value_dict['max_zeta'],value_dict['max_dh'],\
           value_dict['max_eta_pquake'],value_dict['max_speed'],value_dict['max_momentum'],value_dict['max_mflux'],\
           value_dict['tmax'],value_dict['tmin'],value_dict['tfirstpos'],value_dict['tfirstneg'],\
           value_dict['tfirstdraw'],value_dict['tfirstadvis'],value_dict['tfirstwarn']) )


fprint_file.close()


