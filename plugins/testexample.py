import collections
from configuration import hw, sw
import printer
from pprint import pprint
from shuntConvertDatasheet import *

# channel selection
INSPECT_FEDID = 1114
CHANNEL = {"crate":34, "slot":12, "fiber":29, "fibCh":1}
# shunt information
MAX_SHUNT = 31
SHUNT_EVENT_UNIT = 500
SELECTED_SHUNT = [0, 1, 2, 4, 8, 16, 18, 20, 24, 26, 28, 30, 31]
NORM_DIC = {"NORM_FACTR" : 0.0}
ADCtoCharge = shuntConvertDatasheet()

def testexample(raw1={}, raw2={}, book=None, warnQuality=True, fewerHistos=False, **_):
    # sanity check for incoming raw data
    for i, raw in enumerate([raw1, raw2]):
        if not raw:
            continue

        # find the number of time samples in the data
        nTsMax = raw[None]["firstNTs"]
        for fedId, d in sorted(raw.iteritems()):
            if fedId != INSPECT_FEDID:
                continue

            # get the important chunks of the raw data
            blocks = d["htrBlocks"].values()
            header = d["header"]
            evt = header["EvN"]
            shunt = (evt-1)//SHUNT_EVENT_UNIT

            # sanity checks for these chunks
            for block in blocks:
                if type(block) is not dict:
                    printer.warning("FED %d block is not dict" % fedId)
                    continue
                elif "channelData" not in block:
                    printer.warning("FED %d block has no channelData" % fedId)
                    continue
                
                for channelData in block["channelData"].values():  # channel loop
                    # print "evt", evt, block["Slot"], channelData["Fiber"], channelData["FibCh"]
                    # === from here ===
                    if block["Crate"] != 34:
                        continue
                    """
                    fib = 0
                    if block["Slot"] == 12:
                        fib += 12 + channelData["Fiber"] - 1
                        if 13 <= channelData["Fiber"]:  
                            fib -= 2
                    elif block["Slot"] == 11 and 12 <= channelData["Fiber"]:  
                        fib += channelData["Fiber"] - 12 
                    else:  
                        continue """
                    # === to here === > HEP 17 / this result in slot 11 has 0~11 fib and slot 12 has 12~31 fib


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
                            book.fill((ts, adc), "ADC_vs_TS_%s_%d_%d_%d_%d_%d" % (errf, fedId, block["Crate"], block["Slot"], channelData["Fiber"], channelData["FibCh"]),
                                      (nTsMax, nAdcMax), (-0.5, -0.5), (nTsMax - 0.5, nAdcMax - 0.5),
                                      title="FED %d (ErrF %s 0);time slice;ADC;Counts / bin" % (fedId, eq))
                            # fill a 2d histogram where bins on y-axis are ADCs and bins on x-axis are shunt setting 
                            book.fill((shunt, adc), "ADC_vs_Shunt_%s_%d_%d_%d_%d_%d" % (errf, fedId, block["Crate"], block["Slot"], channelData["Fiber"], channelData["FibCh"]),
                                      (nTsMax, nAdcMax), (-0.5, -0.5), (nTsMax - 0.5, nAdcMax - 0.5),
                                      title="FED %d (ErrF %s 0);shunt setting;ADC;Counts / bin" % (fedId, eq))


                        # pedestal subtraction : for each event, the pedestal was estimated by the average of the first two time samples.
                        charge = 0.0
                        pedCharge = 0.0
                        rms_pedCharge = 0.0
                        time_slice = {"ped_ts":[0, 1], "sum_ts":[3,4,5,6]}    # The total charge collected for each channel per event, q, was q = (TS3+TS4+TS5+TS6)-{(TS0+TS1)/2}
                        ped_tsNum = len(time_slice["ped_ts"])
                        sum_tsNum = len(time_slice["sum_ts"])
                        if shunt in SELECTED_SHUNT:
                            for (ts, adc) in enumerate(channelData["QIE"]):
                                if ts in time_slice["ped_ts"]:
                                    pedCharge += float(ADCtoCharge[shunt][adc])
                                    if ts is time_slice["ped_ts"][-1]:
                                        pedCharge /= ped_tsNum
                                if ts in time_slice["sum_ts"]:
                                    charge += float(ADCtoCharge[shunt][adc])
                                    if ts is time_slice["sum_ts"][-1]:
                                        charge -= pedCharge

                            # fill a 2d histogram where bins on y-axis are charges and bins on x-axis are GSel ; pedestal was subtracted 
                            book.fill((shunt, charge), "Charge_vs_GSel_%s_%d_%d_%d_%d_%d" % (errf, fedId, block["Crate"], block["Slot"], channelData["Fiber"], channelData["FibCh"]),
                                      (nTsMax, 200), (-0.5, 0), (nTsMax - 0.5, 130000),
                                      title="FED %d (ErrF %s 0);GSel;Charge[fC];Counts / bin" % (fedId, eq))                                
                            
                            if shunt is 0:
                                # fill mean fedestal
                                book.fill(pedCharge, "pedCharge_%s_%d" % (errf, fedId),
                                          200, 0, 130000,
                                          title="FED %d (ErrF %s 0);pedCharge;Counts / bin" % (fedId, eq))                                




