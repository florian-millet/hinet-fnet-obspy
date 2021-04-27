#!/raid1/fm536/V_envs/Obspy-1-2-2/bin/python3

# packages
import os
import sys
import obspy
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
from   obspy import UTCDateTime
from   obspy.core.inventory import Inventory, Network, Station, Channel, Response

# list of events in GMT time
events = [UTCDateTime(2010,4,11,22,8,11)]

# executable and data folders
datf = '../input_data/fnet/'
sacf = 'sac_files'
ar_t = '.zip'
outf = '../output_data/'

# options
delay   = 0
prepare = True

#------------------#
# prepare the data #
#------------------#
if prepare:

    # loop over files in input folder
    corresp  = [[]]
    zipfiles = sp.Popen('ls '+datf+'*'+ar_t, shell=True, stdout=sp.PIPE)
    for zipout in zipfiles.stdout:
        arch = zipout.decode().strip().split('/')[3]
        print(arch)
    
        # list files inside of archive
        lins = sp.Popen('unzip -l '+datf+arch, shell= True, stdout=sp.PIPE)
        for linz in lins.stdout:
    
            # extract start and end time
            try:
                aa = linz.decode().strip().split()[3].split('_')[2]
                yy = str(aa[ 0: 4]) 
                mo = str(aa[ 4: 6])
                dd = str(aa[ 6: 8])
                hh = str(aa[ 8:10])
                mi = str(aa[10:12])
                ss = str(aa[12:14])
                sta_t = yy+'-'+mo+'-'+dd+'T'+hh+':'+mi+':'+ss
                bb = linz.decode().strip().split()[3].split('_')[3]
                yy = str(bb[ 0: 4]) 
                mo = str(bb[ 4: 6])
                dd = str(bb[ 6: 8])
                hh = str(bb[ 8:10])
                mi = str(bb[10:12])
                ss = str(bb[12:14])
                end_t = yy+'-'+mo+'-'+dd+'T'+hh+':'+mi+':'+ss
                corresp.append([arch,sta_t,end_t])
                break
            except Exception as e:
                pass
    
    # write to file
    c_file = open(datf+'correspondance_file.txt', 'w')
    for element in corresp[1:]:
        for string in element:
            c_file.write(string + ' ')
        c_file.write('\n')
    c_file.close

#------------------#
# loop over events #
#------------------#
for evt in events:

    # create event folder and '.sac' files subfolder in GMT time
    fldr = outf + str(evt)[:19]
    print()
    print(fldr)
    if not os.path.exists(fldr):
        os.makedirs(fldr)
    if not os.path.exists(fldr+'/'+sacf):
        os.makedirs(fldr+'/'+sacf)

    # find correct archive for event
    c_file = open(datf+'correspondance_file.txt', 'r')
    found  = False
    for line in c_file:
        sta_t = UTCDateTime(line.split()[1])
        end_t = UTCDateTime(line.split()[2])
        if (evt+delay >= sta_t) and (evt+delay <= end_t):
            files = line.split()[0]
            found = True
            break
    c_file.close()
    if not found:
        print('Event not found. Try ajusting the \'delay\' parameter')

    # open archive
    sp.call(['unzip', '-qo', datf+files, '-d', fldr+'/'+sacf])

    # loop over channels
    pazfiles = sp.Popen('ls '+fldr+'/'+sacf+'/*.zp', shell=True, stdout=sp.PIPE)
    stream = []
    for pazout in pazfiles.stdout:
        station = pazout.decode().strip().split('/')[4].split('.')[0].split('_')[0]
        channel = pazout.decode().strip().split('/')[4].split('.')[0].split('_')[1]

        # read trace
        trace = obspy.read(fldr+'/'+sacf+'/'+station+'_'+channel+'_*_'+channel)

        # read paz response
        r_file   = open(pazout.decode().strip(),'r')
        pazinfo  = r_file.readline().split()
        constant = float(pazinfo[1])
        pazinfo  = r_file.readline().split()
        n_poles  = int(pazinfo[1])
        poles    = []
        for i in range(n_poles):
            pazinfo  = r_file.readline().split()
            poles.append(complex(float(pazinfo[0]),float(pazinfo[1])))
        pazinfo  = r_file.readline().split()
        n_zeros  = int(pazinfo[1])
        zeros    = []
        for i in range(n_zeros):
            pazinfo  = r_file.readline().split()
            zeros.append(complex(float(pazinfo[0]),float(pazinfo[1])))

        # add response to the data directly
        sts = trace[0].stats
        res = Response.from_paz(zeros=zeros,poles=poles,stage_gain=constant,normalization_frequency=0.02)
        inv = Inventory(networks=[])
        net = Network(code=sts.network,stations=[])
        sta = Station(code=sts.station,latitude=sts.sac['stla'],longitude=sts.sac['stlo'],elevation=sts.sac['stel'])
        cha = Channel(code=sts.channel,latitude=sts.sac['stla'],longitude=sts.sac['stlo'],elevation=sts.sac['stel'],
                      location_code=sts.location,depth=0)
        cha.response=res
        sta.channels.append(cha)
        net.stations.append(sta)
        inv.networks.append(net)
        try:
            trace[0].attach_response(inv)
        except Exception as e:
            print('      THERE WAS AN ERROR WITH THIS TRACE')
            print(e)

        # add data to stream
        if stream == []:
            stream = trace
        else:
            stream.append(trace[0])

    # write in PICKLE format (-fnet)
    name = fldr + '/' + str(evt)[:19] + '-fnet.PICKLE'
    stream.sort()
    stream.write(name, format='PICKLE')





