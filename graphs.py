import math
import configuration
import printer
import utils
r = utils.ROOT()


def combineBinContentAndError(histo, binToContainCombo, binToBeKilled):
    xflows     = histo.GetBinContent(binToBeKilled)
    xflowError = histo.GetBinError(binToBeKilled)

    if xflows == 0.0:  # ugly
        return

    currentContent = histo.GetBinContent(binToContainCombo)
    currentError   = histo.GetBinError(binToContainCombo)

    histo.SetBinContent(binToBeKilled, 0.0)
    histo.SetBinContent(binToContainCombo, currentContent+xflows)

    histo.SetBinError(binToBeKilled, 0.0)
    histo.SetBinError(binToContainCombo, math.sqrt(xflowError**2+currentError**2))


def shiftFlows(histo=None):
    bins = histo.GetNbinsX()
    entries = histo.GetEntries()
    combineBinContentAndError(histo, binToContainCombo=1   , binToBeKilled=0     )
    combineBinContentAndError(histo, binToContainCombo=bins, binToBeKilled=bins+1)
    histo.SetEntries(entries)


def labelXAxis(h=None, labels={}):
    h.SetStats(False)
    axis = h.GetXaxis()
    for iBin in range(1, 1 + axis.GetNbins()):
        x = int(axis.GetBinCenter(iBin))
        axis.SetBinLabel(iBin, labels.get(x, "%d" % x))
    axis.SetLabelSize(1.5*axis.GetLabelSize())
    axis.SetLabelOffset(1.5*axis.GetLabelOffset())


def labelYAxis(h=None, labels={}):
    h.SetStats(False)
    yaxis = h.GetYaxis()
    for iBin, label in sorted(labels.iteritems()):
        yaxis.SetBinLabel(iBin, label)
    yaxis.SetLabelSize(2.0*yaxis.GetLabelSize())


def xMin_xMax(graph=None):
    """assumes that graph is sorted"""
    N = graph.GetN()
    X = graph.GetX()
    if 1 <= N:
        return X[0], X[N-1]
    else:
        printer.error("graph contains zero points.")
        return 0.0, 0.0


def magnify(h, factor=1.0):
    for axis in [h.GetXaxis(), h.GetYaxis()]:
        axis.SetLabelSize(factor*axis.GetLabelSize())
        axis.SetTitleSize(factor*axis.GetTitleSize())


def adjustPad(pad=r.gPad, logY=False, margin=0.2):
    r.gPad.SetLeftMargin(margin)
    r.gPad.SetBottomMargin(margin)
    r.gPad.SetTickx()
    r.gPad.SetTicky()
    if logY:
        r.gPad.SetLogy()


def stylize(h, color=r.kBlack, style=1, width=1):
    #r.gStyle.SetOptStat(110010)
    h.SetStats(False)
    h.SetMinimum(0.5)
    magnify(h, factor=2.0)
    h.SetLineColor(color)
    h.SetLineStyle(style)
    h.SetLineWidth(width)


def legends(legEntries=[]):
    out = []
    dx = 0.8/max(1, len(legEntries))
    x0 = 0.1
    for iLeg, (h, desc) in enumerate(legEntries):
        leg = r.TLegend(x0, 0.91, x0+dx, 1.0)
        x0 += dx
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(h, desc, "l")
        leg.Draw()
        out.append(leg)
    return out


def histoLoop(f, lst, func):
    out = []
    maxes = []
    legEntries = []
    h0 = None
    for i, (x, color, style) in enumerate(lst):
        h = f.Get(func(x))
        if not h:
            continue

        if h.GetName().startswith("ChannelFlavor"):
            labelXAxis(h, labels=configuration.flavorLabels())

        gopts = "hist"
        if i:
            gopts += "same"
        else:
            h0 = h
        maxes.append(h.GetMaximum())
        shiftFlows(h)
        h.Draw(gopts)
        stylize(h, color, style)
        out.append(h)

        legEntries.append((h, h.GetTitle()))
        h.SetTitle("")

    if maxes and h0:
        h0.SetMaximum(2.0*max(maxes))
    out += legends(legEntries)
    return out


def onePage(f=None, pad0=None, pad1=None, pad2=None, feds=[]):
    keep = []

    #label
    pad0.cd()
    text = r.TText(0.5, 0.5, f.GetPath())
    text.SetNDC()
    text.SetTextAlign(22)
    text.SetTextSize(20.0*text.GetTextSize())
    text.Draw()
    keep.append(text)

    #category graphs
    graph = f.Get("category_vs_time")
    if graph:
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(r.gStyle.GetHistLineColor())
        graph.SetMarkerSize(0.5*graph.GetMarkerSize())
        t = graph.GetTitle().split("_")

        pad1.cd()
        adjustPad()
        xMin, xMax = xMin_xMax(graph)
        delta = xMax - xMin
        if delta:
            tenPercent = 0.1 * delta
        else:
            tenPercent = 0.1/60.0  # 1/10 second

        null = r.TH2D("null", ";time (minutes)",
                      1, xMin - tenPercent, xMax + tenPercent,
                      3, 0.5, 3.5)
        null.Draw()
        magnify(null, factor=3.0)
        labelYAxis(null, labels={1: t[0], 2: t[1], 3: t[2]})
        null.GetXaxis().SetTitleOffset(0.75)
        graph.Draw("psame")
        keep += [graph, null]


    #EvN, OrN, BcN agreement
    pad2.cd(1)
    adjustPad(logY=True)

    keep += histoLoop(f,
                      [("OrN", r.kBlue, 1),
                       ("EvN", r.kCyan, 2),
                       ("BcN", r.kBlack, 3),
                       ],
                      lambda x: "delta%s" % x,
                      )

    if True:
        for iHisto, name in enumerate(["", "nBytesSW", "nQieSamples", "nWord16Skipped",
                                       "ErrF1", "ErrF3", "EvN_HTRs", "ChannelFlavor",
                                       "MatchedFibersCh0", "MatchedFibersCh1", "MatchedFibersCh2", "BcN",
                                       #"TTS", "PopCapFrac"
                                       ]):
            if not name:
                continue
            pad2.cd(1+iHisto)
            adjustPad(logY=True)
            h = f.Get(name)
            if h:
                shiftFlows(h)
                h.Draw("hist")
                stylize(h)
                keep.append(h)
            else:
                lst = []
                color = [r.kBlack, r.kRed, r.kBlue]
                color += [r.kBlack] * (len(feds) - len(color))
                style = [1, 2, 3]
                style += [1] * (len(feds) - len(style))
                for iFed, fed in enumerate(sorted(feds)):
                    if name == "BcN" and iFed:
                        continue

                    lst.append((fed, color[iFed], style[iFed]))

                keep += histoLoop(f, lst, lambda x: "%s_%d" % (name, x))

    return keep


def makeSummaryPdf(inputFiles=[], feds=[], pdf="summary.pdf"):
    r.gROOT.SetStyle("Plain")

    canvas = r.TCanvas()
    canvas.Print(pdf + "[")

    pad0 = r.TPad("pad0", "pad0", 0.0, 0.95, 1.0, 1.00)
    pad1 = r.TPad("pad1", "pad1", 0.0, 0.75, 1.0, 0.95)
    pad2 = r.TPad("pad1", "pad1", 0.0, 0.00, 1.0, 0.75)

    pad2.Divide(4, 3)
    pad0.Draw()
    pad1.Draw()
    pad2.Draw()

    for fileName in inputFiles:
        f = r.TFile(fileName)
        if (not f) or f.IsZombie():
            continue
        junk = onePage(f, pad0, pad1, pad2, feds)
        canvas.Print(pdf)
        f.Close()
    canvas.Print(pdf + "]")
