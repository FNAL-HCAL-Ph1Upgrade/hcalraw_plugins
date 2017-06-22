import collections
from configuration import hw, sw
import printer
import math
from pprint import pprint
from shuntConvertDatasheet import *

# channel selection
INSPECT_FEDID = 1114
CHANNEL = {"crate":34, "slot":12, "fiber":14, "fibCh":3}
# shunt information
MAX_SHUNT = 32
SHUNT_EVENT_UNIT = 500
SELECTED_SHUNT = [0, 1, 2, 4, 8, 16, 18, 20, 24, 26, 28, 30, 31]
NORM_DIC = {"NORM_FACTR" : 0.0}
ADCtoCharge = shuntConvertDatasheet()
# pedestal information
PED = { }

def histo_dy(raw1={}, raw2={}, book=None, warnQuality=True, fewerHistos=False, **_):
    # sanity check for incoming raw data
    for i, raw in enumerate([raw1, raw2]):
        if not raw:
            continue

        # find the number of time samples in the data
        nTsMax = raw[None]["firstNTs"]
        for fedId, d in sorted(raw.iteritems()): # event loop
            if fedId != INSPECT_FEDID:
                continue

            # get the important chunks of the raw data
            blocks = d["htrBlocks"].values()
            header = d["header"]
            evt = header["EvN"]
            shunt = (evt-1)//SHUNT_EVENT_UNIT
            print "EVENT", evt
            # sanity checks for these chunks
            for block in blocks:
                if type(block) is not dict:
                    printer.warning("FED %d block is not dict" % fedId)
                    continue
                elif "channelData" not in block:
                    printer.warning("FED %d block has no channelData" % fedId)
                    continue
                # pick a channel >> L48
                # if block["Crate"] != CHANNEL["crate"] or block["Slot"] != CHANNEL["slot"]:
                #    continue
            
                for channelData in block["channelData"].values():  # channel loop
                    # pick a channel >> L43 
                    # if channelData["Fiber"] != CHANNEL["fiber"] or channelData["FibCh"] != CHANNEL["fibCh"]:
                    #    continue
                    
                    # channel information in string form
                    str_ch = "%d_%d_%d_%d"%(block["Crate"], block["Slot"], channelData["Fiber"], channelData["FibCh"])
                    
                    if channelData["QIE"]:
                        # check the error flags
                        errf = "ErrFNZ" if channelData["ErrF"] else "ErrF0"

                        # some logic for naming two histograms, one for a clean error flag and another with a problematic error flag
                        eq = "!=" if channelData["ErrF"] else "=="
                        
                        nAdcMax = 256
                        for (ts, adc) in enumerate(channelData["QIE"]): # TS loop : each channel is consist of 10TS
                            if nTsMax <= ts:
                                break
                            # fill a 2d histogram where bins on y-axis are ADCs and bins on x-axis are time slice index 
                            book.fill((ts, adc), "ADC_vs_TS_%s_%d_%s" % (errf, fedId, str_ch),
                                      (nTsMax, nAdcMax), (-0.5, -0.5), (nTsMax - 0.5, nAdcMax - 0.5),
                                      title="FED %d (ErrF %s 0);time slice;ADC;Counts / bin" % (fedId, eq))
                            if ts is 3: # use max time slice
                                # fill a 2d histogram where bins on y-axis are ADCs and bins on x-axis are shunt setting 
                                book.fill((shunt, adc), "ADC_vs_Shunt_%s_%d_%s" % (errf, fedId, str_ch),
                                          (MAX_SHUNT, nAdcMax), (0, -0.5), (MAX_SHUNT, nAdcMax - 0.5),
                                          title="FED %d (ErrF %s 0);shunt setting;ADC;Counts / bin" % (fedId, eq))


                        # pedestal subtraction : for each event, the pedestal was estimated by the average of the first two time samples.
                        charge = 0.0
                        ADC_pedsub = 0.0
                        pedCharge = 0.0
                        pedADC = 0.0
                        time_slice = {"ped_ts":[1], "sum_ts":[3]}#,4,5,6]}    # The total charge collected for each channel per event, q, was q = (TS3+TS4+TS5+TS6)-{(TS0+TS1)/2}
                        ped_tsNum = len(time_slice["ped_ts"])
                        sum_tsNum = len(time_slice["sum_ts"])
                        if shunt in SELECTED_SHUNT and errf == "ErrF0":
                            for (ts, adc) in enumerate(channelData["QIE"]):
                                if ts in time_slice["ped_ts"]:
                                    pedADC += adc
                                    pedCharge += float(ADCtoCharge[shunt][adc])
                                    if ts is time_slice["ped_ts"][-1]:
                                        pedADC /= ped_tsNum
                                        pedCharge /= ped_tsNum
                                if ts in time_slice["sum_ts"]:
                                    ADC_pedsub += adc
                                    charge += float(ADCtoCharge[shunt][adc])
                                    if ts is time_slice["sum_ts"][-1]:
                                        ADC_pedsub -= pedADC
                                        charge -= pedCharge
                            # accumulate mean value of pedestal for each channel, ( ADC(or charge)sum of corresponding TS / #TS ) 
                            if not(str_ch in PED):
                                PED[str_ch] = []
                            PED[str_ch].append(pedADC) 

                            # fill a 2d histogram where bins on y-axis are charges and bins on x-axis are GSel ; pedestal was subtracted 
                            book.fill((shunt, ADC_pedsub), "ADC_vs_GSel_%s_%d_%s" % (errf, fedId, str_ch),
                                      (MAX_SHUNT, nAdcMax), (0, -0.5), (MAX_SHUNT, nAdcMax-0.5),
                                      title="FED %d (ErrF %s 0);GSel;ADC;Counts / bin" % (fedId, eq))                                
                            book.fill((shunt, charge), "Charge_vs_GSel_%s_%d_%s" % (errf, fedId, str_ch),
                                      (MAX_SHUNT, 200), (0, 20000), (MAX_SHUNT, 100000),
                                      title="FED %d (ErrF %s 0);GSel;Charge[fC];Counts / bin" % (fedId, eq))
        # pedestal mean and RMS 
        if evt == SHUNT_EVENT_UNIT:
            for ch, adclist in PED.items():
                ped_mean = sum(adclist)/len(adclist)
                vsum = 0
                for x in adclist:
                    vsum = vsum + (x-ped_mean)**2
                var = vsum/len(adclist)
                ped_rms = math.sqrt(var)

                book.fill(ped_mean, "mean_pedADC_%s_%d" % (errf, fedId),
                          50, 0, 10,
                          title="FED %d (ErrF %s 0);pedestal mean;Counts / bin" % (fedId, eq)) 
                book.fill(ped_rms, "RMS_pedADC_%s_%d" % (errf, fedId),
                              50, 0, 5,
                              title="FED %d (ErrF %s 0);pedestal RMS;Counts / bin" % (fedId, eq))                             



    
