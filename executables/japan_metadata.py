#!/raid1/fm536/V_envs/Obspy-1-2-2/bin/python3

#----------#
# Preamble #
#----------#
print()

# System integration
import os
import sys
import glob
#import subprocess as sp
import csv

# The usual suspects
import numpy as np
import matplotlib.pyplot as plt

# Obspy - general
import obspy
from   obspy import UTCDateTime

# Obspy - instrument response
from   obspy.core.inventory import Inventory, Network, Station, Channel, Site, Response
from   obspy.io.xseed       import Parser

# Definitions
compute_responses = False #True  #
create_xml        = True  #False #
ev_time = UTCDateTime(2010,4,11,22,8,11) #"2010-04-11T22:08:11")
ev_japn = ev_time + 9*60*60
dat = '../output_data/'
fld = '../responses_as_of_2021-04-19/'
inf = 'dataless_seed/'
ouf = 'resp_files/'
print()
print(str(ev_time)[:19])
print(str(ev_japn)[:19])
print()
if not os.path.exists(fld+ouf):
    os.makedirs(fld+ouf)


#---------------------------------------------------#
# Turn XSEED into STATIONXML for all response files #
#---------------------------------------------------#
if compute_responses:

    # Start with F-NET

    # station info (F-NET)
    with open(fld+inf+'FNET_stations.txt', 'r') as files:
        read = csv.reader(files)
        stal = list(read)
        stal.pop(0)
        for i in range(len(stal)):
            for j in range(len(stal[i])):
                stal[i][j] = stal[i][j].strip()

    # loop over resp files
    for resps in glob.glob(fld+inf+'RESP.B*.txt'): #'RESP.BO.TKA*.txt'): #
        stion = resps.split(sep='/')[3].split(sep='.')[2]
        outpt = fld+ouf+resps.split(sep='/')[3][:-3]+'resp'
        par = Parser(resps)
        blk = par.blockettes

        # number of segments
        for key in blk:
            if key == 50:
                blk_nb = 0
                for kez in blk[key]:
                    blk_nb += 1

        # find station info            
        for i in range(len(stal)):
            if stal[i][1] == stion[:-1]:
                stna = stal[i][0]
                stla = float(stal[i][2][:-1])
                stlo = float(stal[i][3][:-1])
                stel = float(stal[i][4])
                stld = float(stal[i][5]) if stal[i][5] != '' else 0
                print('   ',resps,stion,stna,stla,stlo,stel,stld)
                break

        # rework
        for i in range(blk_nb):
            blk[50][i].network_identifier_code = 1        # Doesn't seem to matter...
            blk[50][i].site_name               = stna     # From channel file
            blk[52][i].instrument_identifier   = 1        # Doesn't seem to matter...
            blk[52][i].latitude                = stla     # From channel file
            blk[52][i].longitude               = stlo     # From channel file
            blk[52][i].elevation               = stel     # From channel file
            blk[52][i].local_depth             = stld     # From channel file
        par.write_xseed(outpt)

    # Then move on to HI-NET

    # loop over resp files
    for resps in glob.glob(fld+inf+'RESP.T*.txt'):
        stion = resps.split(sep='/')[3].split(sep='.')[2]
        outpt = fld+ouf+resps.split(sep='/')[3][:-3]+'resp'
        print('   ', resps)
        par = Parser(resps)
        blk = par.blockettes

        # number of segments
        for key in blk:
            if key == 50:
                blk_nb = 0
                for kez in blk[key]:
                    blk_nb += 1

        # find station info            
        ##for i in range(len(sta2)):
        ##    if sta2[i][1] == stion[:-1]:
        stna = 'type'+stion #sta2[i][0]
        stla = 0.0 #float(sta2[i][2][:-1])
        stlo = 0.0 #float(sta2[i][3][:-1])
        stel = 0.0 #float(sta2[i][4])
        stld = 0.0 #float(sta2[i][5])
        ##        break

        # rework
        # positional info seems not to be used in remov_response, so I hope this is fine...
        for i in range(blk_nb):
            blk[50][i].network_identifier_code = 1        # Doesn't seem to matter...
            blk[50][i].site_name               = stna     # From channel file
            blk[52][i].instrument_identifier   = 1        # Doesn't seem to matter...
            blk[52][i].latitude                = stla     # From channel file
            blk[52][i].longitude               = stlo     # From channel file
            blk[52][i].elevation               = stel     # From channel file
            blk[52][i].local_depth             = stld     # From channel file
        par.write_xseed(outpt)


#-------------------------------------------------#
# Create IRIS-like XML file for the current event #
#-------------------------------------------------#
if create_xml:

    # Load data
    fil_name = dat + str(ev_japn)[:19] + '/' + str(ev_japn)[:19] + '-jpro.PICKLE'
    st = obspy.read(fil_name, format='PICKLE')
    #print(st.__str__(extended=True))
    
    # List networks
    net_all = []
    net_lst = []
    sy = st.copy()
    sy = sy.select(component='Z')
    for tr in st:
        net_all.append(tr.stats.network)
    for elem in net_all:
        if elem not in net_lst:
            net_lst.append(elem)
    
    # List stations within networks
    sta_lst = []
    for elem in net_lst:
        sy = st.copy()
        sy = sy.select(network=elem,component='Z')
        sta_all = []
        for tr in sy:
            sta_all.append(tr.stats.station)
        sta_lst.append(sta_all)
    print()
    print(net_lst)
    print(sta_lst)
    print()
    
    # Create/append xml from scratch
    inv = Inventory(networks = [],
                    source   = "Florian-Millet")
    
    # Loop of networks
    for i in range(len(net_lst)):
        net = Network(code        = net_lst[i],
                      stations    = []) #,
                      #description = "")
    
        # Loop over stations 
        for j in range(len(sta_lst[i])):
            print(sta_lst[i][j], '('+str(sta_lst[i][j][-1:])+')', end=' ... ')
            prb = False
            tra = sy.select(network   = net_lst[i]   ,
                            station   = sta_lst[i][j],
                            component = 'Z')[0]
            sta = Station(code          = tra.stats.station,
                          latitude      = tra.stats.stla   ,
                          longitude     = tra.stats.stlo   ,
                          elevation     = 0.0, #tra.stats.stel   ,
                          creation_date = obspy.UTCDateTime(1980, 1, 1)) #,
                          #site          = Site(name="test station"))
            nwk = 'BO'              if tra.stats.location == '0103' else 'TYPE'
            stn = tra.stats.station if tra.stats.location == '0103' else tra.stats.type
            if stn == 'undefined':
                print('not within the "3 types" advertized, will have to remove later...')
                continue
            res = fld+ouf+'RESP.'+nwk+'.'+stn+'.resp'
            rsp = obspy.read_inventory(res)
            if tra.stats.location == '0103':
                rsp = rsp.select(station=sta_lst[i][j])
    
            # Loop over components
            for k in range(3):
                com = ['Z', 'N', 'E']
                azi = [ 0 ,  0 , 90 ]
                dip = [-90,  0 ,  0 ]
                cha = Channel(code          = "BH"+com[k]   ,
                              location_code = tra.stats.location,
                              latitude      = tra.stats.stla,
                              longitude     = tra.stats.stlo,
                              elevation     = 0.0, #tra.stats.stel,
                              depth         = 0.0           ,
                              azimuth       = azi[k]        ,
                              dip           = dip[k]        ,
                              sample_rate   = tr.stats.sampling_rate)
    
                # Find and append correct response
                if tra.stats.location == '0103':
                    for l in range(len(rsp[0])):
                        chn = rsp[0][l].channels[0]
                        end = UTCDateTime(2599,12,31) if chn.end_date == None else chn.end_date
                        if chn.code == cha.code:
                            if (chn.start_date < ev_japn) and (end > ev_japn):
                                #print('   ',l, chn.code, chn.start_date, end)
                                break
                else:
                    chn = rsp[0][0][0]
                    end = UTCDateTime(2599,12,31) if chn.end_date == None else chn.end_date
                    #print('   ','X', chn.code, chn.start_date, end)
                cha.response = chn.response
    
                # If there is an error for a given station it should be catched here...
                try:
                    sta.channels.append(cha)
                except Exception as e:
                    prb = True

            # And printed here!
            if prb:
                print('had an issue with the response!')
            else:
                print('all went well')

            net.stations.append(sta)
        inv.networks.append(net)
    
    # Finally write to file
    inv_name = dat + str(ev_time)[:19] + '-jap.xml'
    str_name = dat + str(ev_time)[:19] + '-jap.PICKLE'
    inv.write(inv_name, format="stationxml", validate=True)
    st.write( str_name, format="PICKLE")
    


# End
print()
print()




















#__________________________________________________________
# This shouldn't be used but is here for me as a reference
#‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾



## ## ## # Prepare a subset of the data
## ## ## dat_fol = 'OUTPUT_5F_5Hi/' #'output_data/'
## ## ## dat_fmt = '-jap_sta.PICKLE'
## ## ## stream  = obspy.read(dat_fol+str(ev_japn)[:19]+dat_fmt, format='PICKLE')
## ## ## st      = stream.copy()
## ## ## print(st.__str__(extended=True))
## ## ## st      = st.select(station='A*')
## ## ## st.write('RESPONSES_test/2010-jap_data-10_stations.PICKLE', format='PICKLE')


    # THIS FIRST PART IS NOT NEEDED IF I UNDERSTOOD CORRECTLY
    ## # station info (HI-NET)
    ## cntr = 0
    ## temp = []
    ## sta2 = []
    ## with open(fld+'HINET_stations_subset.txt', 'r') as files:
    ##     lines = files.readlines()
    ##     for line in lines:
    ##         cntr += 1
    ##         temp.append(line.split())
    ##         if cntr    ==  3 :
    ##             stna = temp[0][18]
    ##             stco = temp[0][3].split(sep='.')[1][:-1]
    ##             stla = temp[0][13]
    ##             stlo = temp[0][14]
    ##             stel = temp[0][15]
    ##             stld = 0.0
    ##             stty = '1'
    ##             print(temp[0][12])
    ##             if temp[0][12] != '1.023e-07':
    ##                 stty = '2' if temp[0][12] == '1.000e-07' else '3'
    ##             sta2.append([stna,stco,stla,stlo,stel,'','',stty,'',''])
    ##         if line[0] == '#':
    ##             cntr = 0
    ##             temp = []
    ## print(*sta2,sep='\n')


#print()
#print()
#print(inv)
#for i in range(len(inv)):
#    print()
#    print(inv[i])
#    for j in range(len(inv[i])):
#        print()
#        print(inv[i][j])
#        for k in range(len(inv[i][j])):
#            print()
#            print('       |||',i,'|||       |||',j,'|||       |||',k,'|||')
#            print(inv[i][j][k])
#            print(inv[i][j][k].response)



