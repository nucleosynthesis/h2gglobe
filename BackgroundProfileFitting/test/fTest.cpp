#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMsgService.h"
#include "RooDataHist.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = program_options;

bool BLIND = true;

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else {
    cerr << "ERROR -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, string name, double *prob=NULL){
  
  TCanvas *canv = new TCanvas();
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",100,110);
  mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(80),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(80),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(80),Invisible());
  }
  else data->plotOn(plot,Binning(80));

 // data->plotOn(plot,Binning(80));
  pdf->plotOn(plot);
  plot->Draw();
  canv->SaveAs(name.c_str());
  
  if (prob){
	int np = pdf->getParameters(*data)->getSize();
	double chi2 = plot->chiSquare(np);
	*prob = TMath::Prob(chi2*(80-np),80-np);
  }
  
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooAbsData *data, string name, int cat){
  
  int color[5] = {kBlue,kRed,kMagenta,kGreen+1,kYellow+2};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",100,110);
  mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(80),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(80),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(80),Invisible());
  }
  else data->plotOn(plot,Binning(80));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - cat%d",cat),"LEP");
  int i=0;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=5) col=color[i];
    else col=kBlack;
    it->second->plotOn(plot,LineColor(col));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    leg->AddEntry(pdfLeg,it->first.c_str(),"L");
    i++;
  }
  plot->SetTitle(Form("Category %d",cat));
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(name.c_str());
  delete canv;
}

int main(int argc, char* argv[]){
 
  string fileName;
  int ncats;
  string datfile;
  string outDir;
  string outfilename;
  bool is2011=false;
  bool verbose=false;
  bool saveMultiPdf=false;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename)->default_value("multipdfws.root"),         "Save a MultiPdf model with the appropriate pdfs")
    ("is2011",                                                                                  "Run 2011 config")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
  if (vm.count("saveMultiPdf")) saveMultiPdf=true;
  if (vm.count("verbose")) verbose=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }

  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
	outputfile = new TFile(outfilename.c_str(),"RECREATE");
	outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  
  vector<string> functionClasses;
  //functionClasses.push_back("Chebychev");
  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Exponential");
  functionClasses.push_back("PowerLaw");
  functionClasses.push_back("Laurent");
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("Bernstein","pol"));
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("PowerLaw","pow"));
  namingMap.insert(pair<string,string>("Laurent","lau"));
  vector<pair<string,int> > fabChoice;
  if (is2011) {
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
  }
  else {
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
  }

  // store results here
  FILE *resFile = fopen("fTestResults.txt","w");
  vector<map<string,int> > choices_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  pdfsModel.setObsVar(mass);
  
  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
  fprintf(resFile,"\\hline\n");
  for (int cat=0; cat<ncats; cat++){
    
    map<string,int> choices;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    RooDataSet *data = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));

    RooArgList storedPdfs("store");

    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
    fprintf(resFile,"\\hline\n");
    for (vector<string>::iterator funcType=functionClasses.begin(); funcType!=functionClasses.end(); funcType++){
      
      double thisNll=0.;
      double prevNll=0.;
      double chi2=0.;
      double prob=0.;
      int order=1;
      int prev_order=0;
      int cache_order=0;
      RooAbsPdf *prev_pdf;
      RooAbsPdf *cache_pdf;
      while (prob<0.05){
        RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("%d",cat));
        if (!bkgPdf){
          // assume this order is not allowed
          order++;
        }
        else {
          RooFitResult *fitRes = bkgPdf->fitTo(*data,Save(true));
          thisNll = fitRes->minNll();
          //RooAbsReal *nll = bkgPdf->createNLL(*data);
          //RooMinuit m(*nll);
          //m.migrad();
          //thisNll = nll->getVal();
	  double gofProb =0; 
          plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat),&gofProb);
	
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
          prob = TMath::Prob(chi2,order-prev_order);
          cout << "\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          //fprintf(resFile,"%15s && %d && %10.2f && %10.2f && %10.2f \\\\\n",funcType->c_str(),order,thisNll,chi2,prob);
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;

	  if (gofProb > 0.01 && prob<0.05 && saveMultiPdf) {
		std::cout << "Adding " << bkgPdf->GetName() << " "<< gofProb <<std::endl;
      		allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
		storedPdfs.add(*bkgPdf); 
	  }

          prev_order=order;
          prev_pdf=bkgPdf;
          order++;
        }
      }
      fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      choices.insert(pair<string,int>(*funcType,cache_order));
      pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
    }
    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    pdfs_vec.push_back(pdfs);
    plot(mass,pdfs,data,Form("%s/choice_cat%d.pdf",outDir.c_str(),cat),cat);
    plot(mass,pdfs,data,Form("%s/choice_cat%d.png",outDir.c_str(),cat),cat);
    plot(mass,allPdfs,data,Form("%s/multipdf_cat%d.pdf",outDir.c_str(),cat),cat);
    plot(mass,allPdfs,data,Form("%s/multipdf_cat%d.png",outDir.c_str(),cat),cat);

    if (saveMultiPdf){
	  // Put selectedModels into a MultiPdf
	  RooCategory catIndex(Form("pdfindex_%d",cat),"c");
	  RooMultiPdf *pdf = new RooMultiPdf(Form("hgg_cat%d_bkgshape",cat),"all pdfs",catIndex,storedPdfs);
	  RooRealVar nBackground(Form("hgg_cat%d_bkgshape_norm",cat),"nbkg",data->sumEntries(),0,10E8);
	  std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
	  std::cout << "Created MultiPdf in Category " << cat << " total of " << catIndex.numTypes() << std::endl;
	  storedPdfs.Print();
	  std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;

	  mass->setBins(320);
	  RooDataHist dataBinned(Form("binned_data_obs_cat%d",cat),"data",*mass,*data);
	  // Save it (also a binned version of the dataset
	  outputws->import(*pdf);
	  outputws->import(nBackground);
	  outputws->import(catIndex);
	  outputws->import(dataBinned);
    }
    
  }

  FILE *dfile = fopen(datfile.c_str(),"w");
  cout << "Recommended options" << endl;
  for (int cat=0; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",cat); 
    for (map<string,int>::iterator it=choices_vec[cat].begin(); it!=choices_vec[cat].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    fprintf(dfile,"fabian=%s:%d:%s%d\n",fabChoice[cat].first.c_str(),fabChoice[cat].second,namingMap[fabChoice[cat].first].c_str(),fabChoice[cat].second);
    for (map<string,int>::iterator it=choices_vec[cat].begin(); it!=choices_vec[cat].end(); it++){
      if (it->first=="Bernstein"){
        int loword = TMath::Min(3,it->second);
        for (int ord=loword; ord<=fabChoice[cat].second; ord++){
          fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),ord,namingMap[it->first].c_str(),ord);
        }
      }
      else {
        for (int ord=1; ord<=it->second; ord++){
          if ((it->first=="Exponential" || it->first=="PowerLaw") && ord%2==0) continue;
          fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),ord,namingMap[it->first].c_str(),ord);
        }
      }
    }
    fprintf(dfile,"\n");
  }

  inFile->Close();

  if (saveMultiPdf){
	outputfile->cd();
	outputws->Write();
	outputfile->Close();	
  }
  return 0;
}
