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
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <iomanip>
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

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, string name, double *prob){
  
  TCanvas *canv = new TCanvas();
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(80));
  pdf->plotOn(plot_chi2);
  int np = pdf->getParameters(*data)->getSize();
  double chi2 = plot_chi2->chiSquare(np);
  *prob = TMath::Prob(chi2*(80-np),80-np);
  
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
  pdf->paramOn(plot);
  plot->SetTitle("");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("#chi^{2} = %.2f, Prob = %.2f ",chi2*(80-np),*prob));
  canv->SaveAs(name.c_str());
  
  delete canv;
  delete lat;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooAbsData *data, string name, int cat){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
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
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    it->second->plotOn(plot,LineColor(col),LineStyle(style));
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

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}
int getBestFitFunction(RooMultiPdf *bkg, RooAbsData *data, RooCategory *cat, bool silent=false){

	RooRealVar nBackground("bkgshape_norm","nbkg",data->sumEntries(),0,10E8);
	RooExtendPdf extPdf("internal","testextend",*bkg,nBackground);

	double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters(*data);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		std::cout << "CLEAN SET OF PARAMETERS" << std::endl;
		params->Print("V");
		std::cout << "-----------------------" << std::endl;
	}
	
	//bkg->setDirtyInhibit(1);
	RooAbsReal *nllm = extPdf.createNLL(*data,RooFit::Extended());
	RooMinimizer minim(*nllm);
	minim.setStrategy(1);
	
	for (int id=0;id<number_of_indeces;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		if (!silent) {
			std::cout << "BEFORE FITTING" << std::endl;
			params->Print("V");
			std::cout << "-----------------------" << std::endl;
		}
		
		minim.minimize("Minuit2","minimize");
		double minNll = (nllm->getVal())+bkg->getCorrection();
		if (!silent) {
			std::cout << "After Minimization ------------------  " <<std::endl;
			std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
			bkg->Print("v");
			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
			std::cout << " ------------------------------------  " << std::endl;
	
			std::cout << "AFTER FITTING" << std::endl;
			params->Print("V");
			std::cout << " Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
			std::cout << " Correction Applied is " << bkg->getCorrection() <<std::endl;
			std::cout << " NLL + c = " <<std::setprecision(10) <<  minNll << std::endl;
			std::cout << "-----------------------" << std::endl;
		}
			
		if (minNll < global_minNll){
        		global_minNll = minNll;
			snap.assignValueOnly(*params);
        		best_index=id;
		}
	}
	params->assignValueOnly(snap);
    	cat->setIndex(best_index);
	
	if (!silent) {
		std::cout << "Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
		bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
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
		("unblind", 																																								"Dont blind plots")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
	if (vm.count("unblind")) BLIND=false;
  if (vm.count("saveMultiPdf")) {
	saveMultiPdf=true;
  } else {
	saveMultiPdf=false;
  }
  if (vm.count("verbose")) verbose=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }

  std::cout << "SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
	outputfile = new TFile(outfilename.c_str(),"RECREATE");
	outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");

	transferMacros(inFile,outputfile);
	RooRealVar *intL = (RooRealVar*)inWS->var("IntLumi");
	outputws->import(*intL);
  
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
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
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
		// now need more randoms as have 15 cats
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
  }

  // store results here
  FILE *resFile = fopen("fTestResults.txt","w");
  vector<map<string,int> > choices_vec;
  vector<map<string,std::vector<int> > > choices_envelope_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  pdfsModel.setObsVar(mass);
  double upperEnvThreshold = 0.1; // upper threshold on delta(chi2) to include function in envelope (looser than truth function)
  
  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
  fprintf(resFile,"\\hline\n");
  for (int cat=0; cat<ncats; cat++){
    
    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    RooDataSet *dataFull = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));
    mass->setBins(320);

    RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
    RooDataSet *data = (RooDataSet*)&thisdataBinned;
   

    RooArgList storedPdfs("store");

    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
    fprintf(resFile,"\\hline\n");

    double MinimimNLLSoFar=1e10;
    int simplebestFitPdfIndex = 0;
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
      std::vector<int> pdforders;
      while (prob<0.05){
        RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order);
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
	  double fail=0;
          plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat),&fail);
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
          prob = TMath::Prob(chi2,order-prev_order);
          cout << "\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          //fprintf(resFile,"%15s && %d && %10.2f && %10.2f && %10.2f \\\\\n",funcType->c_str(),order,thisNll,chi2,prob);
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;
          prev_order=order;
          prev_pdf=bkgPdf;
          order++;
        }
      }
      fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      choices.insert(pair<string,int>(*funcType,cache_order));
      pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

      // Now run loop to determine functions inside envelope
      if (true/*saveMultiPdf*/){
        chi2=0.;
        thisNll=0.;
        prevNll=0.;
        prob=0.;
        order=1;
        prev_order=0;
        cache_order=0;
	std::cout << "Determining Envelope Functions for Family " << *funcType << ", cat " << cat << std::endl;
	std::cout << "Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;
        while (prob<upperEnvThreshold){
         RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("pdf_%d",cat));
          if (!bkgPdf ){
          // assume this order is not allowed
          order++;
          }

          else {
           RooFitResult *fitRes = bkgPdf->fitTo(*data,Save(true));
           thisNll = fitRes->minNll();
	   double myNll = 2.*thisNll;
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

	   if (gofProb > 0.01 && prob < upperEnvThreshold ) { // Looser requirements for the envelope
		std::cout << "Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
			  << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
      		allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
		storedPdfs.add(*bkgPdf);
		pdforders.push_back(order);
 
		if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
		 // Pval based correction
		 simplebestFitPdfIndex = storedPdfs.getSize()-1;
		 MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
		}

	   }

           prev_order=order;
           prev_pdf=bkgPdf;
           order++;
        }
      }
      
      fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    }

    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);
    plot(mass,pdfs,data,Form("%s/truths_cat%d.pdf",outDir.c_str(),cat),cat);
    plot(mass,pdfs,data,Form("%s/truths_cat%d.png",outDir.c_str(),cat),cat);
    plot(mass,allPdfs,data,Form("%s/multipdf_cat%d.pdf",outDir.c_str(),cat),cat);
    plot(mass,allPdfs,data,Form("%s/multipdf_cat%d.png",outDir.c_str(),cat),cat);

    if (saveMultiPdf){
	  // Put selectedModels into a MultiPdf
	  std::string ext = is2011 ? "7TeV" : "8TeV";
	  RooCategory catIndex(Form("pdfindex_%d_%s",cat,ext.c_str()),"c");
	  RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_cat%d_%s_bkgshape",cat,ext.c_str()),"all pdfs",catIndex,storedPdfs);
	
	  RooRealVar nBackground(Form("CMS_hgg_cat%d_%s_bkgshape_norm",cat,ext.c_str()),"nbkg",data->sumEntries(),0,10E8);
	  //double check the best pdf!
	  int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,true);
	  catIndex.setIndex(bestFitPdfIndex);
	  std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
	  std::cout << "Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
	  storedPdfs.Print();
	  std::cout << "Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
	  std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
	  std::cout << " Simple check of index "<< simplebestFitPdfIndex <<std::endl;

	  mass->setBins(320);
	  //RooDataHist dataBinned(Form("binned_data_obs_cat%d",cat),"data",*mass,*data);
	  RooDataHist dataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);

	  // Save it (also a binned version of the dataset
	  outputws->import(*pdf);
	  outputws->import(nBackground);
	  outputws->import(catIndex);
	  outputws->import(dataBinned);
		outputws->import(*data);
    }
    
  }
  if (saveMultiPdf){
	outputfile->cd();
	outputws->Write();
	outputfile->Close();	
  }

  //FILE *dfile = fopen(datfile.c_str(),"w");
  cout << "Recommended options" << endl;
  /*
  for (int cat=0; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",cat); 
    for (map<string,int>::iterator it=choices_vec[cat].begin(); it!=choices_vec[cat].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    fprintf(dfile,"fabian=%s:%d:%s%d\n",fabChoice[cat].first.c_str(),fabChoice[cat].second,namingMap[fabChoice[cat].first].c_str(),fabChoice[cat].second);
    for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat].begin(); it!=choices_envelope_vec[cat].end(); it++){
	std::vector<int> ords = it->second;
        for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
         // if ((it->first=="Exponential" || it->first=="PowerLaw") && ord%2==0) continue;
          fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
        }
     // }
    }
    fprintf(dfile,"\n");
  }
  */

  inFile->Close();

  return 0;
}
