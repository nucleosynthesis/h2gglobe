#!/usr/bin/env python

import os
import sys
import shlex

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datfile",dest="datfile",help="Read from datfile")
parser.add_option("-f","--file",dest="files",default=[],action="append",help="Add a file")
parser.add_option("-o","--outname",dest="outname",help="Name of output pdf/png/C")
parser.add_option("-c","--color",dest="colors",default=[],action="append",help="Set color")
parser.add_option("-s","--style",dest="styles",default=[],action="append",help="Set style")
parser.add_option("-w","--width",dest="widths",default=[],action="append",help="Set width")
parser.add_option("-n","--name",dest="names",default=[],action="append",help="Set name")
parser.add_option("-t","--text",dest="text",type="string",default="",help="Add Text")
parser.add_option("-e","--expected",dest="expected",default=False,action="store_true",help="Expected only")
parser.add_option("-m","--method",dest="method",type="string",help="Method to run")
parser.add_option("-l","--legend",dest="legend",type="string",help="Legend position - x1;y1;x2;y2")
parser.add_option("-y","--yaxis",dest="yaxis",type="string",help="y-axis range - y1;y2")
parser.add_option("","--limit",dest="limit",default=False,action="store_true",help="Do limit plot")
parser.add_option("","--pval",dest="pval",default=False,action="store_true",help="Do p-value plot")
parser.add_option("","--maxlh",dest="maxlh",default=False,action="store_true",help="Do best fit mu plot")
parser.add_option("","--mh",dest="mh",default=False,action="store_true",help="Do NLL mass scan plot")
parser.add_option("","--mu",dest="mu",default=False,action="store_true",help="Do NLL mu scan plot")
parser.add_option("","--mumh",dest="mumh",default=False,action="store_true",help="Do NLL mu vs mh scan plot")
parser.add_option("","--rvrf",dest="rvrf",default=False,action="store_true",help="Do NLL rv vs rf scan plot")
parser.add_option("-v","--verbose",dest="verbose",default=False,action="store_true")
parser.add_option("-b","--batch",dest="batch",default=False,action="store_true")
(options,args)=parser.parse_args()

# Required for back compatbility:
if options.limit: options.method='limit'
if options.pval: options.method='pval'
if options.maxlh: options.method='maxlh'
if options.mh: options.method='mh'
if options.mu: options.method='mu'
if options.mumh: options.method='mumh'
if options.rvrf: options.method='rvrf'
if not options.outname: options.outname=options.method

allowed_methods=['pval','limit','maxlh','mh','mu','mumh','rvrf']
if not options.datfile and options.method not in allowed_methods:
  print 'Invalid method. Must set one of: ', allowed_methods
  sys.exit()

import ROOT as r
outf = r.TFile('CombinePlotCanvases.root','RECREATE')

if options.batch: r.gROOT.SetBatch()

# make canvas and dummy hist
canv = r.TCanvas()
canv.SetGrid(True)
dummyHist = r.TH1D("dummy","",1,110,150)
dummyHist.GetXaxis().SetTitle('m_{H} (GeV)')
dummyHist.SetTitleSize(.05,"X")
dummyHist.SetTitleOffset(0.75,"X")
dummyHist.SetTitleSize(.05,"Y")
dummyHist.SetTitleOffset(0.75,"Y")
# make text box
lat = r.TLatex()
lat.SetNDC()
lat.SetTextSize(0.03)
lat.SetTextFont(42);

def pvalPlot(allVals):
  
  canv.Clear()
  if options.verbose: print 'Plotting pvalue...'
  canv.SetLogy(True)
  mg = r.TMultiGraph()
  if not options.legend: leg = r.TLegend(0.6,0.12,0.89,0.4)
  else: leg = r.TLegend(float(options.legend.split(';')[0]),float(options.legend.split(';')[1]),float(options.legend.split(';')[2]),float(options.legend.split(';')[3]))
  leg.SetFillColor(0)

  # make graphs from values
  for k, values in enumerate(allVals):
    graph = r.TGraph()
    for j in range(len(values)):
      graph.SetPoint(j,values[j][0],values[j][1])
      print "Add point ", values[j][0], values[j][1]
    
    graph.SetLineColor(int(options.colors[k]))
    graph.SetLineStyle(int(options.styles[k]))
    graph.SetLineWidth(int(options.widths[k]))
    if options.names[k]!="-1": leg.AddEntry(graph,options.names[k],'L')
    mg.Add(graph)
 
  # draw dummy hist and multigraph
  dummyHist.GetYaxis().SetTitle('Local p-value')
  mg.Draw("A")
  if not options.yaxis:
    dummyHist.SetMinimum(mg.GetYaxis().GetXmin())
    dummyHist.SetMaximum(mg.GetYaxis().GetXmax())
  else:
    dummyHist.SetMinimum(float(options.yaxis.split(';')[0]))
    dummyHist.SetMaximum(float(options.yaxis.split(';')[1]))
    
  dummyHist.SetLineColor(0)
  dummyHist.SetStats(0)
  dummyHist.Draw("AXIS")
  mg.Draw("L")
  dummyHist.Draw("AXIGSAME")

  # draw sigma lines
  sigmas=[1,2,3,4,5]
  lines=[]
  labels=[]
  for i,sig in enumerate(sigmas):
    y = r.RooStats.SignificanceToPValue(sig)
    if options.verbose: print sig, y
    lines.append(r.TLine(110,y,150,y))
    lines[i].SetLineWidth(2)
    lines[i].SetLineStyle(2)
    lines[i].SetLineColor(r.kRed)
    labels.append(r.TLatex(110 + 2, y * 1.1, "%d #sigma" % (i+1)))
    labels[i].SetTextAlign(11);
    if y<=dummyHist.GetMaximum() and y>=dummyHist.GetMinimum():
      lines[i].Draw('SAME')
      labels[i].Draw('SAME')

  # draw text
  lat.DrawLatex(0.12,0.92,"CMS Preliminary")
  lat.DrawLatex(0.75,0.94,options.text)
  
  # draw legend
  leg.Draw()
  canv.RedrawAxis()

  # print canvas
  canv.Update()
  if not options.batch: raw_input("Looks ok?")
  canv.Print('%s.pdf'%options.outname)
  canv.Print('%s.png'%options.outname)
  canv.Print('%s.C'%options.outname)
  canv.SetName(options.outname)
  outf.cd()
  canv.Write()

def maxlhPlot(allVals):

  canv.Clear()
  canv.SetLogy(False)
  if options.verbose: print 'Plotting maxlh...'
  mg = r.TMultiGraph()
  if not options.legend: leg = r.TLegend(0.6,0.8,0.89,0.89)
  else: leg = r.TLegend(float(options.legend.split(';')[0]),float(options.legend.split(';')[1]),float(options.legend.split(';')[2]),float(options.legend.split(';')[3]))
  leg.SetFillColor(0)
 
  # make graph from values
  for k, values in enumerate(allVals):
    graph = r.TGraphAsymmErrors()
    point_counter=0
    for j in range(len(values)):
      if (j%4==0):
        mh = values[j][0]
        fit = values[j][1]
        low = values[j+1][1]
        high = values[j+2][1]
        graph.SetPoint(point_counter,mh,fit)
        graph.SetPointError(point_counter,0,0,abs(fit-low),abs(high-fit))
        point_counter+=1
        if options.verbose: print mh, fit, low, high
    
    graph.SetMarkerStyle(21)
    graph.SetMarkerSize(0.5)
    graph.SetLineColor(1)
    graph.SetLineWidth(2)
    graph.SetFillColor(r.kGreen)
    leg.AddEntry(graph,'68% CL Band','F')
    mg.Add(graph)
  
  # draw dummy hist and multigraph
  dummyHist.GetYaxis().SetTitle('Best fit #sigma/#sigma_{SM}')
  mg.Draw("A")
  if not options.yaxis:
    dummyHist.SetMinimum(mg.GetYaxis().GetXmin())
    dummyHist.SetMaximum(mg.GetYaxis().GetXmax())
  else:
    dummyHist.SetMinimum(float(options.yaxis.split(';')[0]))
    dummyHist.SetMaximum(float(options.yaxis.split(';')[1]))
  dummyHist.SetLineColor(0)
  dummyHist.SetStats(0)
  dummyHist.Draw("AXIS")
  mg.Draw("3")
  mg.Draw("LPX")
  dummyHist.Draw("AXIGSAME")
  
  # draw line at y=1 
  l = r.TLine(110,1.,150,1.)
  l.SetLineColor(r.kBlue)
  l.SetLineStyle(r.kDashed)
  l.SetLineWidth(2)
  l.Draw()
  
  # draw line at y=0 
  l2 = r.TLine(110,0.,150,0.)
  l2.SetLineColor(r.kRed)
  l2.SetLineStyle(r.kDashed)
  l2.SetLineWidth(2)
  l2.Draw()
  
  # draw text
  lat.DrawLatex(0.12,0.92,"CMS Preliminary")
  lat.DrawLatex(0.75,0.94,options.text)
  
  # draw legend
  leg.Draw()
  canv.RedrawAxis()

  # print canvas
  canv.Update()
  if not options.batch: raw_input("Looks ok?")
  canv.Print('%s.pdf'%options.outname)
  canv.Print('%s.png'%options.outname)
  canv.Print('%s.C'%options.outname)
  canv.SetName(options.outname)
  outf.cd()
  canv.Write()

def limitPlot(allVals):

  canv.Clear()
  canv.SetLogy(False)
  if options.verbose: print 'Plotting limit...'
  mg = r.TMultiGraph()
  if not options.legend: leg = r.TLegend(0.6,0.7,0.89,0.89)
  else: leg = r.TLegend(float(options.legend.split(';')[0]),float(options.legend.split(';')[1]),float(options.legend.split(';')[2]),float(options.legend.split(';')[3]))
  leg.SetFillColor(0)

  # make graph from values
  for k, values in enumerate(allVals):
    graph = r.TGraphAsymmErrors()
    exp = r.TGraphAsymmErrors()
    oneSigma = r.TGraphAsymmErrors()
    twoSigma = r.TGraphAsymmErrors()
    point_counter=0
    for j in range(len(values)):
      if (j%6==0):
        mh = values[j][0]
        down95 = values[j][1]
        down68 = values[j+1][1]
        median = values[j+2][1]
        up68 = values[j+3][1]
        up95 = values[j+4][1]
        obs = values[j+5][1]
        graph.SetPoint(point_counter,mh,obs)
        exp.SetPoint(point_counter,mh,median)
        oneSigma.SetPoint(point_counter,mh,median)
        oneSigma.SetPointError(point_counter,0,0,abs(median-down68),abs(up68-median))
        twoSigma.SetPoint(point_counter,mh,median)
        twoSigma.SetPointError(point_counter,0,0,abs(median-down95),abs(up95-median))
        point_counter+=1
        if options.verbose: print mh, median, down68, up68, down95, up95, obs
    
    graph.SetMarkerStyle(21)
    graph.SetMarkerSize(0.5)
    graph.SetLineColor(1)
    graph.SetLineWidth(2)
    exp.SetLineColor(1)
    exp.SetLineStyle(2)
    oneSigma.SetLineStyle(2)
    twoSigma.SetLineStyle(2)
    oneSigma.SetFillColor(r.kGreen)
    twoSigma.SetFillColor(r.kYellow)
    if len(allVals)>1:
      exp.SetLineColor(int(options.colors[k]))
      exp.SetLineStyle(2)
      exp.SetLineWidth(int(options.widths[k]))
      graph.SetMarkerColor(int(options.colors[k]))
      graph.SetLineColor(int(options.colors[k]))
      leg.AddEntry(graph,options.names[k],'L')
    else:
      exp.SetLineColor(1)
      exp.SetLineStyle(2)
      exp.SetLineWidth(2)
      if not options.expected: leg.AddEntry(graph,'Observed','L')
      #leg.AddEntry(exp,'Expected','L')
      leg.AddEntry(oneSigma,'Expected #pm 1#sigma','FL') 
      leg.AddEntry(twoSigma,'Expected #pm 2#sigma','FL') 
    
    if len(allVals)==1:
      mg.Add(twoSigma)
      mg.Add(oneSigma)
    mg.Add(exp)
    if not options.expected: mg.Add(graph)
  
  # draw dummy hist and multigraph
  dummyHist.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%%CL} / \sigma(H#rightarrow #gamma #gamma)_{SM}")
  mg.Draw("A")
  if not options.yaxis:
    dummyHist.SetMinimum(mg.GetYaxis().GetXmin())
    dummyHist.SetMaximum(mg.GetYaxis().GetXmax())
  else:
    dummyHist.SetMinimum(float(options.yaxis.split(';')[0]))
    dummyHist.SetMaximum(float(options.yaxis.split(';')[1]))
  dummyHist.SetLineColor(0)
  dummyHist.SetStats(0)
  dummyHist.Draw("AXIS")
  mg.Draw("3")
  mg.Draw("LPX")
  dummyHist.Draw("AXIGSAME")
 
  # draw line at y=1 
  l = r.TLine(110,1.,150,1.)
  l.SetLineColor(r.kRed)
  l.SetLineWidth(2)
  l.Draw()

  # draw text
  lat.DrawLatex(0.12,0.92,"CMS Preliminary")
  lat.DrawLatex(0.75,0.94,options.text)
  
  # draw legend
  leg.Draw()
  canv.RedrawAxis()

  # print canvas
  canv.Update()
  if not options.batch: raw_input("Looks ok?")
  canv.Print('%s.pdf'%options.outname)
  canv.Print('%s.png'%options.outname)
  canv.Print('%s.C'%options.outname)
  canv.SetName(options.outname)
  outf.cd()
  canv.Write()

def runStandard():
  config = []
  for k, f in enumerate(options.files):
    tf = r.TFile(f)
    tree = tf.Get('limit')
    values=[]
    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      values.append([tree.mh, tree.limit])
    values.sort(key=lambda x: x[0])
    config.append(values)

  if options.method=='pval': pvalPlot(config)
  elif options.method=='limit': limitPlot(config)
  elif options.method=='maxlh': maxlhPlot(config)

def read1D(file,x,i,xtitle):
  tree = file.Get('limit')
  tree.Draw('2*deltaNLL:%s'%x,'','')
  gr = r.gROOT.FindObject('Graph').Clone('gr_%s_%d'%(x,i))
  gr.SetTitle("")
  gr.GetXaxis().SetTitle(xtitle)
  gr.GetYaxis().SetTitle("-2 #Delta LL")

  gr.Sort()
  last = None
  for i in range(gr.GetN(),0,-1):
    if gr.GetX()[i-1] == last:
      gr.RemovePoint(i-1)
    last = gr.GetX()[i-1]
  return gr

def plot1DNLL():
 
  canv.Clear()
  canv.SetLogy(False)
  if options.method=='mh':
    x = 'MH'
    xtitle = 'm_{H} (GeV)'
  elif options.method=='mu':
    x = 'r'
    xtitle = '#sigma / #sigma_{SM}'

  canv = r.TCanvas(x,x,500,500)
  if not options.legend: leg  = r.TLegend(0.35,0.65,0.6,0.79)
  else: leg = r.TLegend(float(options.legend.split(';')[0]),float(options.legend.split(';')[1]),float(options.legend.split(';')[2]),float(options.legend.split(';')[3]))
  leg.SetLineColor(0)
  leg.SetFillColor(0)

  mg = r.TMultiGraph()
  objs = []
  graphs = []
  clean_graphs = []
  for k, f in enumerate(options.files):
    tf = r.TFile(f)
    graph = read1D(tf,x,k,xtitle)
    graph.SetLineColor(int(options.colors[k]))
    graph.SetLineStyle(int(options.styles[k]))
    graph.SetLineWidth(int(options.widths[k]))
    clean_graph = graph.Clone()
    clean_graphs.append(clean_graph)
    leg.AddEntry(clean_graph,options.names[k],'L')

    # convert graph into function
    sp = r.GraphToTF1('mygraph%d'%k,graph)
    func = r.TF1("myfunc%d" %k,sp,graph.GetX()[0],graph.GetX()[graph.GetN()-1],1,"GraphToTF1")
    func.SetParameter(0,0.)
    func.SetLineColor(0)
    func.SetLineStyle(0)
    graphs.append(graph)
    objs.append(sp)
    objs.append(func)
    objs.append(graph)
  
    xmin = func.GetMinimumX()
    eminus = xmin - func.GetX(1.,func.GetXmin(),xmin)
    eplus  = func.GetX(1.,xmin,func.GetXmax()) - xmin
    eminus2 = xmin - func.GetX(4.,func.GetXmin(),xmin)
    eplus2  = func.GetX(4.,xmin,func.GetXmax()) - xmin
    print "%s : %1.4f +%1.3g -%1.3g" % ( graph.GetName(), xmin, eplus , eminus )

    # for the first passed file only get the intersection lines
    if k==0:
      rng = 6
      axmin = 999.
      axmax = -999.
      fit = xmin
      err = (abs(eplus)+abs(eminus))/2.
      eminusR = xmin - func.GetX(rng,func.GetXmin(),xmin)
      eplusR  = func.GetX(rng,xmin,func.GetXmax()) - xmin
      axmin = min(axmin,xmin - eminusR)
      axmax = max(axmax,xmin + eplusR)
      lines = [ r.TLine(axmin, 0, axmax, 0),
                r.TLine(axmin, 1, axmax, 1), r.TLine(xmin-eminus,  0, xmin-eminus,  1), r.TLine(xmin+eplus,  0, xmin+eplus,  1), 
                r.TLine(axmin, 4, axmax, 4), r.TLine(xmin-eminus2, 0, xmin-eminus2, 4), r.TLine(xmin+eplus2, 0, xmin+eplus2, 4) ]
    
  dH = r.TH1D("dH","",1,axmin,axmax)
  dH.GetXaxis().SetTitle(xtitle)
  dH.SetTitleSize(.04,"X")
  dH.SetTitleOffset(0.95,"X")
  dH.SetTitleSize(.04,"Y")
  dH.SetTitleOffset(0.85,"Y")
  if options.method=='mh': dH.GetXaxis().SetNdivisions(505)
  dH.GetYaxis().SetTitle('-2 #Delta LL')
  if not options.yaxis: dH.GetYaxis().SetRangeUser(0.,rng)
  else: dH.GetYaxis().SetRangeUser(float(options.yaxis.split(';')[0]),float(options.yaxis.split(';')[1]))
  dH.SetLineColor(0)
  dH.SetStats(0)
  dH.Draw("AXIS")
    
  for gr in clean_graphs:
    gr.GetXaxis().SetRangeUser(axmin,axmax)
    gr.GetXaxis().SetNdivisions(505)
    if not options.yaxis: gr.GetYaxis().SetRangeUser(0.,rng)
    else: gr.GetYaxis().SetRangeUser(float(options.yaxis.split(';')[0]),float(options.yaxis.split(';')[1]))
    gr.Draw("L")

  # draw legend
  if len(options.files)>1:
    leg.Draw('same')
        
  # draw intersection lines
  for l in lines:
    l.SetLineWidth(2)
    l.SetLineColor(r.kRed)
    l.Draw('same')
  
  # draw text
  lat.DrawLatex(0.12,0.92,"CMS Preliminary")
  lat.DrawLatex(0.7,0.94,options.text)

  # draw fit value
  lat2 = r.TLatex()
  lat2.SetNDC()
  lat2.SetTextSize(0.04)
  lat2.SetTextAlign(22)
  if options.method=='mh': lat2.DrawLatex(0.5,0.85,"m_{H} = %5.1f #pm %3.1f"%(fit,err))
  elif options.method=='mu': lat2.DrawLatex(0.5,0.85,"#sigma/#sigma_{SM} = %3.1f #pm %3.1f"%(fit,err))

  canv.Update()
  if not options.batch: raw_input("Looks ok?")
  canv.Print('%s.pdf'%options.outname)
  canv.Print('%s.png'%options.outname)
  canv.Print('%s.C'%options.outname)
  canv.SetName(options.outname)
  outf.cd()
  canv.Write()

def plot2DNLL(col,type,xtitle,ytitle):
  canv.Clear()
  canv.SetLogy(False)
  canv = r.TCanvas("%s"%type,"%s"%type,500,500)
  tf = r.TFile('%s_bands.root'%type)
  if not options.legend: leg = r.TLegend(0.7,0.7,0.88,0.88)
  else: leg = r.TLegend(float(options.legend.split(';')[0]),float(options.legend.split(';')[1]),float(options.legend.split(';')[2]),float(options.legend.split(';')[3]))
  leg.SetFillColor(0)
  
  gbest = tf.Get('%s_best'%type)
  th2 = tf.Get('%s_th2'%type)
  cont68 = tf.Get('%s_c68'%type)
  cont95 = tf.Get('%s_c95'%type)
  r.gStyle.SetOptStat(0)
  th2.SetTitle("")
  th2.GetXaxis().SetTitle(xtitle)
  th2.GetYaxis().SetTitle(ytitle)
  th2.GetYaxis().SetTitleOffset(1.2)
  th2.GetXaxis().SetTitleSize(0.04)
  th2.GetYaxis().SetTitleSize(0.04)
  if options.method=='mumh':
    th2.GetXaxis().SetRangeUser(122,128)
    if not options.yaxis: th2.GetYaxis().SetRangeUser(0,2.5)
    else: th2.GetYaxis().SetRangeUser(float(options.yaxis.split(';')[0]),float(options.yaxis.split(';')[1]))
  if options.method=='rvrf':
    th2.GetXaxis().SetRangeUser(-1.5,4.5)
    if not options.yaxis: th2.GetYaxis().SetRangeUser(-1.5,4.5)
    else: th2.GetYaxis().SetRangeUser(float(options.yaxis.split(';')[0]),float(options.yaxis.split(';')[1]))
  th2.GetZaxis().SetRangeUser(0,10)
  if col: th2.Draw("colz")
  else: th2.Draw("axis")


  gbest.Draw("Psame")
  for i in range(cont68.GetSize()):
    g68 = cont68.At(i).Clone("")
    g68.SetLineWidth(2)
    g68.SetLineColor(1)
    g68.Draw("csame")
  for i in range(cont95.GetSize()):
    g95 = cont95.At(i).Clone("")
    g95.SetLineWidth(2)
    g95.SetLineColor(1)
    g95.SetLineStyle(2)
    g95.Draw("csame")
 
  leg.AddEntry(gbest,"Best Fit","PL");
  leg.AddEntry(g68,"1#sigma","L");
  leg.AddEntry(g95,"2#sigma","L");
  leg.Draw("same")

  # draw text
  lat.DrawLatex(0.12,0.92,"CMS Preliminary")
  lat.DrawLatex(0.7,0.94,options.text)

  #write results
  if options.method=='mumh':
    mh=gbest.GetX()[0]
    mu=gbest.GetY()[0]
    latexmu = r.TLatex()
    latexmu.SetNDC()
    latexmu.SetTextSize(0.035)
    latexmu.DrawLatex(0.12,0.84,"#mu = %5.2f"%(mu))
    latexmh = r.TLatex()
    latexmh.SetNDC()
    latexmh.SetTextSize(0.035)
    latexmh.DrawLatex(0.12,0.80,"m_{H} = %5.1f GeV"%(mh))  
  if options.method=='rvrf':
    rf=gbest.GetX()[0]
    rv=gbest.GetY()[0]
    latexrv = r.TLatex()
    latexrv.SetNDC()
    latexrv.SetTextSize(0.035)
    latexrv.DrawLatex(0.67,0.22,"#mu_{qqH+VH} = %5.2f "%(rv))
    latexrf = r.TLatex()
    latexrf.SetNDC()
    latexrf.SetTextSize(0.035)
    latexrf.DrawLatex(0.67,0.18,"#mu_{ggH+ttH} = %5.2f "%(rf)) 
    
  
  if col:
    canv.Update()
    if not options.batch: raw_input("Looks ok?")
    canv.Print('%s_col.pdf'%options.outname)
    canv.Print('%s_col.png'%options.outname)
    canv.Print('%s_col.C'%options.outname)
    canv.SetName('%s_col'%options.outname)
    outf.cd()
    canv.Write()
  else:
    canv.Update()
    if not options.batch: raw_input("Looks ok?")
    canv.Print('%s.pdf'%options.outname)
    canv.Print('%s.png'%options.outname)
    canv.Print('%s.C'%options.outname)
    canv.SetName(options.outname)
    outf.cd()
    canv.Write()

def run2DNLL():
  if len(options.files)!=1: 
    sys.exit('For the 2D Scan plots you should only pass one file!')
  r.gROOT.LoadMacro('makeBands.cxx')
  if options.method=='mumh': 
    tf = r.TFile('MuMH_bands.root','RECREATE')
    r.readMassScan2D(tf,"MuMH",options.files[0])
    tf.Close()
  if options.method=='rvrf': 
    tf = r.TFile('RvRf_bands.root','RECREATE')
    r.readParamScan2D(tf,"RvRf",options.files[0],"RF","RV")
    tf.Close()

  if options.method=='mumh': 
    plot2DNLL(True,'MuMH','m_{H} (GeV)','#sigma / #sigma_{SM}')
    plot2DNLL(False,'MuMH','m_{H} (GeV)','#sigma / #sigma_{SM}')
  if options.method=='rvrf': 
    plot2DNLL(True,'RvRf','#mu_{ggH+ttH}','#mu_{qqH+VH}')
    plot2DNLL(False,'RvRf','#mu_{ggH+ttH}','#mu_{qqH+VH}')

def run():
  if options.verbose:
    print options.method
    print options.files
    print options.colors
    print options.styles
    print options.widths
    print options.names

  if options.method=='pval' or options.method=='limit' or options.method=='maxlh':
    runStandard()
  elif options.method=='mh' or options.method=='mu':
    r.gROOT.ProcessLine(".x rootPalette.C")
    r.gROOT.LoadMacro('GraphToTF1.C+')
    plot1DNLL()
  elif options.method=='mumh' or options.method=='rvrf':
    r.gROOT.ProcessLine(".x rootPalette.C")
    run2DNLL()

# __MAIN__

if options.datfile:
  d = open(options.datfile)
  for line in d.readlines():
    if line.startswith('#'): continue
    if line=='\n': continue
    config={}
    line = line.replace('\=','EQUALS')
    for opt in line.split(':'):
      config[opt.split('=')[0]] = opt.split('=')[1].replace('EQUALS','=').split(',')
    for opt in ['colors','styles','widths']:
      if opt in config.keys():
        config[opt] = [int(x) for x in config[opt]]
    
    for key, item in config.items():
      if len(item)==1 and key in ['method','text','outname','legend','yaxis']:
        item=item[0].strip('\n')
      setattr(options,key,item)

    if options.verbose: print options
    run()

else:
  run()

print 'All canvases written to:', outf.GetName()
outf.Close()
