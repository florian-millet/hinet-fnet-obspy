#!/raid1/fm536/V_envs/Obspy-1-2-2/bin/python3

import os
import sys
import obspy
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
from   obspy import UTCDateTime

# list of events in GMT time
events = [UTCDateTime(2010,4,11,22,8,11)]

# travel time calculations
from obspy.taup import TauPyModel
from obspy      import geodetics
pha = ['ScS','Sdiff']
mod = TauPyModel('prem')
ede = 619.6
ela = 37.0
elo = -3.5

# executable and data folders
exef = '../../../executables/'  #'../../../../win32/win32tools/'
conv = 'w32tow1'                #'w32tow1.src/w32tow1'
datf = '../input_data/'
winf = 'win_files'
ar_t = '.zip' #'.tar.gz'
outf = '../output_data/'

# options
plot_yes_no = False

#------------------#
# loop over events #
#------------------#
for evt in events:

    # create event folder and '.win' files subfolder in JST time
    ev2  = evt  + 9*60*60
    fldr = outf + str(ev2)[:19]
    print()
    print(fldr)
    if not os.path.exists(fldr):
        os.makedirs(fldr)
    if not os.path.exists(fldr+'/'+winf):
        os.makedirs(fldr+'/'+winf)

    # initialize prm files list ~ will first contain all files then no duplicate files
    prm_list = []

    #---------------------------------------------------#
    # loop over X minutes archive files for given event #
    #---------------------------------------------------#
    ev       = str(ev2)[:4] + str(ev2)[5:7] + str(ev2)[8:10]
    tarfiles = sp.Popen('ls '+datf+'*'+ev+'*'+ar_t, shell=True, stdout=sp.PIPE)
    for tarout in tarfiles.stdout:
    
        # create folder for X minutes archive and unzip said archive
        arch   = tarout.decode().strip().split('/')[2]
        code   = arch[0:2]+arch[3:5]
        archiv = arch[6:10]+'-'+arch[10:12]+'-'+arch[12:14]+'T'+arch[14:16]+':'+arch[16:18]
        print(' ', archiv, code, end=' ')
        if not os.path.exists(fldr+'/'+archiv):
            os.makedirs(fldr+'/'+archiv)
        if ar_t == '.tar.gz':
            sp.call(['tar', '-C', fldr+'/'+archiv+'/', '-xzf', datf+arch])
        if ar_t == '.zip':
            sp.call(['unzip', '-qo', datf+arch, '-d', fldr+'/'+archiv+'/'])
        os.chdir(fldr+'/'+archiv)

        # read individual unziped prm files
        prmfiles = sp.Popen('ls *'+ev+'*.ch' , shell=True, stdout=sp.PIPE)
        for prmout in prmfiles.stdout:
            prms = prmout.decode().strip()
            with open(prms, 'r') as pfile:
                linz = pfile.readlines()
                for lin in linz:
                    if lin.split()[0] == '#':
                        continue
                    prm_list.append([lin.split()])

        # convert win32 to win for obspy
        winfiles = sp.Popen('ls '+ev+'*.cnt', shell=True, stdout=sp.PIPE)
        print(' at  ',*os.getcwd().split('/')[8:],sep='/')
        for winout in winfiles.stdout:
            win32 = winout.decode().strip()
            sp.call(['./' + exef + conv, '-s', win32, '../' + winf + '/' + win32[:-4] + '.win' ])

        # safely go to next subfolder
        os.chdir('../../')

    # remove duplicates in prm files
    prm_nodup = []
    for elem in prm_list:
        if elem not in prm_nodup:
            prm_nodup.append(elem)
    prm_list = prm_nodup
    #print(*prm_list,sep='\n')

    #------------------------------------#
    # read and merge resulting win files #
    #------------------------------------#
    infiles = sp.Popen('ls '+fldr+'/'+winf+'/'+ev+'*.win', shell=True, stdout=sp.PIPE)
    st      = []
    for inout in infiles.stdout:

        # read individual files
        current = inout.decode().strip()
        code    = current[-10:-6]
        print('   ', current, code)
        cu = obspy.read(current, format='WIN')

        # make sure that 'dtype' is the same (i.e. convert to float64)
        for tr in cu:
            tr.data = np.float64(tr.data)
            tr.stats.location = code

        # append and merge
        st = cu.copy() if len(st)==0 else st+cu
        try:
            st.merge()
        except Exception as e:
            print(e)
            sys.exit()

    # write to file
    name = fldr + '/' + str(ev2)[:19] + '-jraw.PICKLE'
    st.sort()
    st.write(name, format='PICKLE')

    #------------------------#
    # write stuff into stats #
    #------------------------#
    name = fldr + '/' + str(ev2)[:19] + '-jraw.PICKLE'
    st = obspy.read(name, format='PICKLE')
    for tr in st:

        # find correct metadata parameter line
        found = False
        for prm_elem in prm_list:
            metadata = prm_elem[0]
            if tr.stats.channel == metadata[0]:
                found = True
                break
        if not found:
            print('channel not found')

        # start extracting metadata
        lat = metadata[13]
        lon = metadata[14]
        sen = metadata[7]
        res = metadata[8]
        hxe = metadata[0]
        val = metadata[7]
        try:
            net = metadata[3].split(sep='.')[0]
            sta = metadata[3].split(sep='.')[1]
        except Exception as e:
            net = 'X'
            sta = metadata[3]

        # type of logger for HiNet stations
        typ = 'individual'
        if metadata[3][-1:] == 'H':
            if   metadata[12] == '1.023e-07':
                typ = '1'
            elif metadata[12] == '1.000e-07':
                typ = '2'
            elif metadata[12] == '1.192e-07':
                typ = '3'
            else:
                typ = 'undefined'

        # channel orientation is a mess...
        if len(metadata[4]) > 1:
            ind_channel = metadata[4][0] if metadata[4][1] in ['B'] else metadata[4][1]
        else:
            ind_channel = metadata[4]
        if ind_channel not in ['N', 'E', 'Z', 'X', 'Y', 'U', 'V']:
            st.remove(tr)
            #print('ATTENTION', tr, ind_channel)
        cha = 'BH' + ind_channel
        cha = 'BHZ' if cha in ['BHU','BHV'] else cha
        cha = 'BHN' if cha in ['BHX']       else cha
        cha = 'BHE' if cha in ['BHY']       else cha

        # write metadata
        tr.stats.network     = net
        tr.stats.station     = sta
        tr.stats.channel     = cha
        tr.stats.stla        = lat
        tr.stats.stlo        = lon
        tr.stats.sensitivity = sen
        tr.stats.response    = res 
        #tr.stats.location    = hxe
        tr.stats.type        = typ
        tr.stats.correction  = str(float(val)/200)

    # remove when duplicates or missing components
    st.sort()
    cur_sta = 'ABCDEF'
    cntr_tot = 0
    cntr_pbl = 0
    cntr_dup = 0
    for tr in st:
        if cur_sta != tr.stats.station:
            cntr_tot += 1
            cur_sta = tr.stats.station
            sy = st.select(station=cur_sta)
            if   len(sy)   == 3:
                pass
            elif len(sy)%3 == 0:
                cntr_dup += 1
                sy = st.select(station=cur_sta)
                for rmv in sy[3:]:
                    st.remove(rmv)
                continue
            elif len(sy)   != 3:
                cntr_pbl += 1
                sy = st.select(station=cur_sta)
                for rmv in sy:
                    st.remove(rmv)
                continue
    print('problems   =', cntr_pbl)
    print('duplicates =', cntr_dup)
    print('total      =', cntr_tot)
    print('good       =', len(st)/3, '('+str(len(st))+')')

    # write to file
    name = fldr + '/' + str(ev2)[:19] + '-jpro.PICKLE'
    st.write(name, format='PICKLE')

    #----------------------#
    # plot map of stations #
    #----------------------#

    #--------------------------#
    # plot individual stations #
    #--------------------------#
    if plot_yes_no :
        st = obspy.read(fldr+'/' + str(ev2)[:19]+'-jpro.PICKLE', format='PICKLE')
        st.sort()

        st = st.select(component='Z')
        for tr in st:
            print()
            print()
            print(tr.stats.station[-1:])
            print(tr.stats)
        print()
        print()
        sys.exit()

        startt = st[0].stats.starttime
        st.trim(startt+600-300,startt+600+3300)
        #print(st.__str__(extended=True))
        colors_1 = ['0101','0103','0301','0203']
        colors_2 = [   'r',   'k',   'b',   'g']
        if not os.path.exists('plots/'+fldr):
            os.makedirs('plots/'+fldr)
        sta_list = []
        for tr in st:
            if tr.stats.station not in sta_list:
                sta_list.append(tr.stats.station)
        for sta in sta_list:
            # select and filter
            sy = st.select(station=sta)
            cd = sy[0].stats.location
            co = colors_2[colors_1.index(cd)] if cd in colors_1 else 'c'
            nb = len(sy)
            startt = sy[0].stats.starttime
            sy.filter('bandpass',freqmin=0.03,freqmax=0.09,corners=4,zerophase=True)
            # compute ttimes
            sla = sy[0].stats.stla
            slo = sy[0].stats.stlo
            dis = geodetics.gps2dist_azimuth(ela,elo,float(sla),float(slo))[0]/1000/111
            tim = mod.get_travel_times(source_depth_in_km=ede,
                                       distance_in_degree=dis,
                                       phase_list=pha)
            # rotate (?)
            baz = geodetics.gps2dist_azimuth(ela,elo,float(sla),float(slo))[2]
            sy.rotate('NE->RT',back_azimuth=baz)
            sy.sort()
            # differentiate and flip(?)
            sy.differentiate()
            fl = -1
            # trim
            width = 200
            delay = 300
            times = [tim[0].time+delay-width, tim[0].time+delay+width]
            sz = sy.copy()
            sz.trim(startt+times[0],startt+times[1])
            # plot
            plt.figure(1,figsize=(20,10),dpi=100)
            for i in range(nb):
                plt.subplot(nb,2,2*i+1) #i+1)
                plt.plot(sy[i].times(),fl*sy[i].data,co)
                #plt.plot(sy[i].times('matplotlib'),sy[i].data)
                plt.legend(sy[i].stats.channel[2])
                plt.title(sy[0].stats.station)
                amp = np.max(np.abs(sy[i].data))
                plt.plot([times[0],times[0]],[-amp,amp])
                plt.plot([times[1],times[1]],[-amp,amp])
            for i in range(nb):
                plt.subplot(nb,2,2*i+2)
                plt.plot(sz[i].times(),fl*sz[i].data,co)
                #plt.plot(sz[i].times('matplotlib'),sz[i].data)
                plt.legend(sz[i].stats.channel[2])
                plt.title(sz[0].stats.station)
            print('   ',sy[0].stats.station,'   ',dis,'   ',tim[0])
            fnam = 'plots/'+fldr+'/'+sy[0].stats.station+'-3.png'
            plt.savefig(fnam)
            plt.close()
#-------------------------
            #plt.show()
            #sys.exit()
#-------------------------
        
    




print()
