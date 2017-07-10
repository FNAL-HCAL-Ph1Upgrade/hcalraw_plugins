import collections
from configuration import hw, sw
import printer
import numpy as np
from pprint import pprint
from shuntConvertDatasheet import *

""" =========================================================================
CHECK LIST ( search KEY#)
0. SWITCH : KEY0
1. INSPECT_FEDID : KEY1
2. CHANNEL : KEY2
3. SHUNT_EVENT_UNIT and MAX_SHUNT : KEY3
4. PEDESTAL : KEY4
========================================================================= """
# switch >> KEY0
PLOT_SWITCH = {"ADC_vs_TS":"ON", "Charge_vs_TS":0, "ADC_vs_Shunt":0, "ADC_vs_GSel":0, "Charge_vs_GSel":0, "pedMean":0, "pedRMS":0, "Phase":"ON"}

# channel selection
INSPECT_FEDID = 1114 # >>> KEY1
CHANNEL = {"crate":34, "slot":12, "fiber":20, "fibCh":0} # >>> KEY2
MAX_PHASE = 100
PHASE_EVENT_UNIT = 100
# shunt information
MAX_SHUNT = 32
SHUNT_EVENT_UNIT = 100000 # ex) if it is 500 then, evt1~500 = shunt 0, evt501~1000 = shunt1, etc.. >>> KEY3
SELECTED_SHUNT = [0, 1, 2, 4, 8, 16, 18, 20, 24, 26, 28, 30, 31]
LINEARIZATION = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 11.5]
ADCtoCharge = shuntConvertDatasheet()
# pedestal information 
PED = { }
# The total charge collected for each channel per event, q, was q = (TS3+TS4+TS5+TS6)-{(TS0+TS1)/2} >>> KEY4
time_slice = {"ped_ts":[1], "sum_ts":[2]}#,4,5,6]}    

def makeHisto(raw1={}, raw2={}, book=None, warnQuality=True, fewerHistos=False, **_):
    # sanity check for incoming raw data
    for i, raw in enumerate([raw1, raw2]):
        if not raw:
            continue
        # find the number of time samples in the data
        nTsMax = raw[None]["firstNTs"]
        for fedId, d in sorted(raw.iteritems()): # event loop
            if fedId == None : # >>> KEY1
                continue

            # get the important chunks of the raw data
            blocks = d["htrBlocks"].values()
            header = d["header"]
            evt = header["EvN"]
            shunt = (evt-1)//SHUNT_EVENT_UNIT 
            shunt = SELECTED_SHUNT[shunt] # activate this line when shunt scan runs over only Gsel, KEY3
            phase = (evt-1)//PHASE_EVENT_UNIT
            if evt%SHUNT_EVENT_UNIT == 0: 
                print "EVENT", evt
            # sanity checks for these chunks
            for block in blocks:
                if type(block) is not dict:
                    printer.warning("FED %d block is not dict" % fedId)
                    continue
                elif "channelData" not in block:
                    printer.warning("FED %d block has no channelData" % fedId)
                    continue
                # pick a channel >>> KEY2
                #if block["Crate"] != CHANNEL["crate"] or block["Slot"] != CHANNEL["slot"]:
                #    continue
            
                for channelData in block["channelData"].values():  # channel loop
                    # pick a channel >> KEY2 
                    #if channelData["Fiber"] >= 8 : # != CHANNEL["fiber"] or channelData["FibCh"] != CHANNEL["fibCh"]:
                    #    continue
                    
                    # channel information in string form
                    str_ch = "%d_%d_%d_%d"%(block["Crate"], block["Slot"], channelData["Fiber"], channelData["FibCh"])
                    
                    if channelData["QIE"]:
                        # check the error flags
                        errf = "ErrFNZ" if channelData["ErrF"] else "ErrF0"
                        # some logic for naming two histograms, one for a clean error flag and another with a problematic error flag
                        eq = "!=" if channelData["ErrF"] else "=="

                        # delet this if you want to make plots of channel with error
                        #if errf == "ErrFNZ":
                        #    continue
                         
                        nAdcMax = 256
                        ts_max = 2 # you can find this time slice by seeing pulse shape plot whici is ADC vs TS (or charge vs TS)
                        for (ts, adc) in enumerate(channelData["QIE"]): # TS loop : each channel is consist of 10TS
                            if nTsMax <= ts:
                                break
                            
                            if PLOT_SWITCH["ADC_vs_TS"] is "ON":
                                # fill a 2d histogram where bins on y-axis are ADCs and bins on x-axis are time slice index 
                                book.fill((ts, adc), "ADC_vs_TS_%s_%d_%s" % (errf, fedId, str_ch),
                                          #(nTsMax, nAdcMax), (-0.5, -0.5), (nTsMax - 0.5, nAdcMax - 0.5),    # active this to make 2D histogram
                                          nTsMax, -0.5, nTsMax - 0.5,                                         # active this to make TProfile
                                          title="FED %d_%s (ErrF %s 0);time slice;ADC;Counts / bin" % (fedId, str_ch, eq))
                            if PLOT_SWITCH["Charge_vs_TS"] is "ON":
                                # fill a TProfile where bins on y-axis are Charge
                                if shunt in SELECTED_SHUNT:
                                    charge = float(ADCtoCharge[shunt][adc])
                                    book.fill((ts, charge), "Charge_vs_TS_%s_%d_%s" % (errf, fedId, str_ch),
                                              nTsMax, -0.5, nTsMax-0.5,
                                              title="FED %d_%s (ErrF %s 0);time slice;Charge[fC];Counts / bin" % (fedId, str_ch , eq))
                            if PLOT_SWITCH["ADC_vs_Shunt"] is "ON":
                                if ts is ts_max: # use max time slice
                                    # fill a 2d histogram where bins on y-axis are ADCs and bins on x-axis are shunt setting 
                                    book.fill((shunt, adc), "ADC_vs_Shunt_%s_%d_%s" % (errf, fedId, str_ch),
                                              (MAX_SHUNT, nAdcMax), (0, -0.5), (MAX_SHUNT, nAdcMax - 0.5),
                                              title="FED %d_%s (ErrF %s 0);shunt setting;ADC;Counts / bin" % (fedId, str_ch, eq))
                            if PLOT_SWITCH["Phase"] is "ON" :
                                if ts == 3:
                                    deno = adc
                                    if deno == 0:
                                        print  evt, str_ch
                                if ts == 4 and deno != 0:
                                    nume = adc
                                    book.fill((phase, nume/deno), "TS4TS3_vs_phasesetting_%d_%s" % ( fedId, str_ch),
                                              MAX_PHASE, 0, MAX_PHASE,
                                              title="FED %d_%s ;phase setting;TS4/TS3;Counts / bin" % (fedId, str_ch))

                        charge_pedsub = 0.0
                        ADC_pedsub = 0.0
                        pedCharge = 0.0
                        pedADC = 0.0
                        ped_tsNum = len(time_slice["ped_ts"])
                        sum_tsNum = len(time_slice["sum_ts"])
                        if shunt in SELECTED_SHUNT :#and errf != "ErrF0":
                            # convert shunt to Gsel
                            gsel = SELECTED_SHUNT.index(shunt)
                            #print evt, shunt, LINEARIZATION[gsel]
                            for (ts, adc) in enumerate(channelData["QIE"]):
                                # accumulate mean value of pedestal for each channel, ( ADC(or charge)sum of corresponding TS / #TS ) 
                                if ts == (3):
                                    if not(str_ch in PED):
                                        PED[str_ch] = []
                                    PED[str_ch].append(adc) 
                                # pedestal subtraction
                                if ts in time_slice["ped_ts"]:
                                    pedADC += adc
                                    pedCharge += float(ADCtoCharge[1][adc]) # 1/x : ADCtoCharge[1][adc] otherwise ADCtoCharge[shunt][adc] KEY
                                    if ts is time_slice["ped_ts"][-1]:
                                        pedADC /= ped_tsNum
                                        pedCharge /= ped_tsNum
                                if ts in time_slice["sum_ts"]:
                                    ADC_pedsub += adc
                                    charge_pedsub += float(ADCtoCharge[1][adc]) # 1/x : ADCtoCharge[1][adc] otherwise ADCtoCharge[shunt][adc] KEY
                                    if ts is time_slice["sum_ts"][-1]:
                                        ADC_pedsub -= pedADC
                                        charge_pedsub -= pedCharge
                            
                            # fill 2D histograms where bins on y-axis are charges and bins on x-axis are GSel ; pedestal was subtracted 
                            # if it is shunt scan fill 'ADC_pedsub' and charge_pedsub' / otherwise fill 'pedADC' and 'pedCharge' >>> KEY4 
                            if PLOT_SWITCH["ADC_vs_GSel"] is "ON":
                                book.fill((LINEARIZATION[gsel], ADC_pedsub), "ADC_vs_ShuntValue_%d_%s" % (fedId, str_ch),
                                          # (MAX_SHUNT, nAdcMax), (0, -0.5), (MAX_SHUNT, nAdcMax-0.5),                    # active this to make 2D histogram 
                                          #(24, nAdcMax), (0, -0.5), (12, nAdcMax-0.5),                    # active this to make 2D histogram 
                                          #MAX_SHUNT, 0, MAX_SHUNT,
                                          24, 0, 12,
                                          title="FED %d_%s;shunt value;ADC;Counts / bin" % (fedId, str_ch))                                
                            if PLOT_SWITCH["Charge_vs_GSel"] is "ON":
                                #print "#evt", evt, str_ch, " ",LINEARIZATION[gsel], charge_pedsub
                                book.fill((LINEARIZATION[gsel], charge_pedsub), "Charge_vs_ShuntValue_%d_%s" % (fedId, str_ch),
                                          #(24, 500), (0, 0), (12, 150000),                            # active this to make 2D histogram 
                                          24, 0, 12,
                                          title="FED %d_%s;shunt value;Charge[fC];Counts / bin" % (fedId, str_ch))
                                          #title="FED %d_%s (ErrF %s 0);GSel;events;Counts / bin" % (fedId, str_ch, eq))

            #print "# of channel", len(PED)

            # pedestal mean and RMS 
            if evt == SHUNT_EVENT_UNIT and (PLOT_SWITCH["pedMean"] is "ON" or PLOT_SWITCH["pedRMS"] is "ON"):
                for ch, adclist in PED.items():
                    #print "channel : ", ch, "lenth : ", len(adclist)
                    ped_mean = np.mean(adclist)
                    ped_rms = np.std(adclist)
                    if ped_rms == 0:
                        print "The Signal of channel ", ch," is flat. This channel won't be considered." 
                        continue
                
                    # fill 1D histograms 
                    if PLOT_SWITCH["pedMean"] is "ON":
                        book.fill(ped_mean, "mean_pedADC_%s_%d" % (errf, fedId),
                                  50, 0, 10,
                                  title="FED %d;pedestal mean;Counts / bin" % (fedId)) 
                    if PLOT_SWITCH["pedRMS"] is "ON":
                        book.fill(ped_rms, "RMS_pedADC_%d" % (fedId),
                                  100, 0, 1,
                                  title="FED %d;pedestal RMS;Counts / bin" % (fedId))


    
