#include "PdfWeightSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include "TMath.h"
#include <assert.h>
#include <algorithm>

#include "Macros/Normalization_8TeV.h"

PdfWeightSmearer::PdfWeightSmearer(const  std::string & theFile,  Normalization_8TeV * norm, std::string dId , std::string uId ) : 
	efficiency_file(theFile), norm_(norm), downId(dId), upId(uId)
{
  name_="PdfWeightSmearer";
}

PdfWeightSmearer::~PdfWeightSmearer()
{
}

void PdfWeightSmearer::readFile(std::string uId, std::string dId ){

  kFactorSmearers_.resize(6);

  // first is central/up/down for the G-fusion, then its the VBF
  TH2F* temp = (TH2F*) thePdfWeightFile_->Get("GF_cent");
  assert(temp!=0);    kFactorSmearers_[0]=(TH2F*) temp->Clone(("Hmasscent")); kFactorSmearers_[0]->SetDirectory(0);

  temp = (TH2F*) thePdfWeightFile_->Get( Form("GF_%s",uId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[1]=(TH2F*) temp->Clone(("Hmass_up")); kFactorSmearers_[1]->SetDirectory(0);

  temp = (TH2F*) thePdfWeightFile_->Get( Form("GF_%s",dId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[2]=(TH2F*) temp->Clone(("Hmass_down")); kFactorSmearers_[2]->SetDirectory(0);

  // Need to normalize to central inumber of entries!
  double C_integral = kFactorSmearers_[0]->GetEntries();

  kFactorSmearers_[0]->Scale(C_integral/kFactorSmearers_[0]->Integral());
  kFactorSmearers_[2]->Scale(C_integral/kFactorSmearers_[2]->Integral());
  kFactorSmearers_[1]->Scale(C_integral/kFactorSmearers_[1]->Integral());

  temp = (TH2F*) thePdfWeightFile_->Get("VBF_cent");
  assert(temp!=0);    kFactorSmearers_[3]=(TH2F*) temp->Clone(("Hmasscent_vbf")); kFactorSmearers_[3]->SetDirectory(0);

  temp = (TH2F*) thePdfWeightFile_->Get( Form("VBF_%s",uId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[4]=(TH2F*) temp->Clone(("Hmass_up_vbf")); kFactorSmearers_[4]->SetDirectory(0);

  temp = (TH2F*) thePdfWeightFile_->Get( Form("VBF_%s",dId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[5]=(TH2F*) temp->Clone(("Hmass_down_vbf")); kFactorSmearers_[5]->SetDirectory(0);

  // Need to normalize to central number of entries integral!
  C_integral = kFactorSmearers_[3]->GetEntries();
  kFactorSmearers_[3]->Scale(C_integral/kFactorSmearers_[3]->Integral());
  kFactorSmearers_[5]->Scale(C_integral/kFactorSmearers_[5]->Integral());
  kFactorSmearers_[4]->Scale(C_integral/kFactorSmearers_[4]->Integral());
}

bool PdfWeightSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  if( sample_type >= 0 ) { return true; }
  int genMassPoint = std::round(norm_->GetMass(sample_type));

  int process_shift=0;
  if ( norm_->GetProcess(sample_type) == "ggh" ){
  	process_shift =  0;
  } else if (norm_->GetProcess(sample_type) == "vbf"){
  	process_shift =  3;
  } else {
    return true;
  }

  if( genMassPoint > 150 ) { genMassPoint=150; } // Warning: missing k-factor
  if( genMassPoint == 100 ) { genMassPoint=105; }  // Warning: missing k-factor

 
  double kWeight = getWeight( p4, nPu, syst_shift , process_shift);
  weight = (kWeight > 0) ? kWeight : 0;
  return true;
}


bool PdfWeightSmearer::init() 
{
  cout << name_ << " - Opening PdfWeight file"<<endl;
  thePdfWeightFile_ = TFile::Open( efficiency_file.c_str() );
  assert(thePdfWeightFile_!=0);

  readFile(upId, downId);
  thePdfWeightFile_->Close();

  return true;
}

double PdfWeightSmearer::getPdfWeight(int genMassPoint, int id, double gPT , double gY ) const 
{
    const TH2F* tmp = kFactorSmearers_[id]; 
    return tmp->GetBinContent(tmp->FindFixBin(gY,gPT));
    assert(0); return 0.;
}

/*
double pvalPoisDiff(int k1, int ksum){
	//  Eqn 2.1-2.2 from  http://www.ucs.louisiana.edu/~kxk4695/JSPI-04.pdf
	
	// run the sums and get the min of two 
	int k2 = ksum-k1;
	double pc = 0.5; // n,lamba same for each 
	// since assime rate is same for this hypothesis and MC sample is the same for each 
	double pk = TMath::Binomial(ksum,k1)*TMath::Power(pc,k1)*TMath::Power(1-pc,k2);
	
	double ps = 0.;
	if (k1 > 100) return 1.;

	for (int i=k1+1;i<=ksum;i++){
	  ps+=TMath::Binomial(ksum,i)*TMath::Power(pc,i)*TMath::Power(1-pc,ksum-i);
	} 
	return ps+pk < 1-ps+pk ? ps+pk :1-ps+pk ;
}

bool checkConsistency(double xd, double yd ){

	// assume x/y
	if (y==0) return true;
	int x = (int) xd;
	int y = (int) yd;

	// Check p-val  
	double pval = pvalPoisDiff(x,x+y);
	if (pval > 0.05) return true;
	return false; 
	
} 
*/

double PdfWeightSmearer::getWeight( const TLorentzVector & p4, const int nPu, float syst_shift, int process_shift) const
{
  float gPT = p4.Pt();
  float gY  = fabs(p4.Rapidity());

  if (gPT > 300) return 1.;
  if (gY  > 5) return 1.;

  int    varId     = syst_shift > 0 ?   1 : 2;
  varId+=process_shift; // this is a ggH or qqH signal

  double nominal   = getPdfWeight( 0, process_shift, gPT, gY );
  double variation = getPdfWeight( 0, varId, gPT, gY );


  return ( max( 1. + fabs(syst_shift) * ((variation/nominal)-1), 0.) ) ;

}
