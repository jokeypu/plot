#include <iostream>
using namespace std;

void ReadBump()
{
  TFile *inputFile = new TFile("digi_ScanMom3.root");
  TTree* inputTree = (TTree*)inputFile->Get("cbmsim");
  if (!inputTree) {
    cout<<"Error: cannot find Input TTree"<<endl;
	return;
  }

  TClonesArray* fBumpArray = new TClonesArray("PndEmcBump");
  inputTree->SetBranchAddress("EmcBump", &fBumpArray);
  TClonesArray* fSharedDigiArray = new TClonesArray("PndEmcSharedDigi");
  inputTree->SetBranchAddress("EmcSharedDigi", &fSharedDigiArray);


  vector<TH2D*> h2DBumpPosition;
  int nevents = inputTree->GetEntriesFast();
  nevents = 1;
  // loop events
  for (int ievt = 0; ievt<nevents; ievt++){
    inputTree->GetEntry(ievt);
	int nbumps = fBumpArray->GetEntriesFast();
	cout<<"The number of bump is "<<nbumps<<endl;

    // loop bumps
	for (int ibump = 0; ibump<nbumps; ibump++) {
	  PndEmcBump* theBump = (PndEmcBump*)fBumpArray->At(ibump);
	  std::vector<Int_t> digiList = theBump->DigiList();

	  string hname = Form("h2_%2d", ibump);
	  h2DBumpPosition.push_back(new TH2D(hname.c_str(), "", 71, 0, 71, 160, 0, 160));
	  int ih2 = h2DBumpPosition.size()-1;
	  int digiSize = digiList.size();
	  // create 2D histogram for each bump;
	  for (int idigi=0; idigi<digiSize; idigi++) {
	    PndEmcDigi* digi = (PndEmcDigi*)fSharedDigiArray->At(idigi);
		double energy = digi->GetEnergy();
		int iTheta = digi->GetThetaInt();
		int iPhi   = digi->GetPhiInt();
		h2DBumpPosition.at(ih2)->Fill(iTheta, iPhi, 10+log(energy));
	    //cout<<"log energy "<<log(energy)<<endl;
	  }
	}
  }

  TCanvas *c1 = new TCanvas("c1","");
  c1->SetMargin(0.15, 0.05, 0.15, 0.05);
  gStyle->SetOptStat(0);
  int nh2 = h2DBumpPosition.size();
  if (nh2<1) return ;
  for (int ih2=0; ih2<nh2; ih2++) {
    h2DBumpPosition.at(ih2)->SetFillColor(ih2+1);
    if (ih2==0) {
	  h2DBumpPosition.at(ih2)->GetXaxis()->SetTitle("Theta ID");
	  h2DBumpPosition.at(ih2)->GetYaxis()->SetTitle("Phi ID");
	  h2DBumpPosition.at(ih2)->GetXaxis()->CenterTitle();
	  h2DBumpPosition.at(ih2)->GetYaxis()->CenterTitle();
	  h2DBumpPosition.at(ih2)->GetXaxis()->SetTitleSize(0.06);
	  h2DBumpPosition.at(ih2)->GetYaxis()->SetTitleSize(0.06);
	  h2DBumpPosition.at(ih2)->GetXaxis()->SetLabelSize(0.06);
	  h2DBumpPosition.at(ih2)->GetYaxis()->SetLabelSize(0.06);
	  h2DBumpPosition.at(ih2)->Draw("box");
    }
	else h2DBumpPosition.at(ih2)->Draw("boxsame");
  }
  return ;
}
