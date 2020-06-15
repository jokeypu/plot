#include "bes3plotstyle.C"
void plot6(){
	TH1D* h1 = new TH1D("h1", "e1", 100, 0.0, 0.25);
	TH1D* h2 = new TH1D("h2", "e9", 100, 0.0, 0.25);
	TH1D* h3 = new TH1D("h3", "e25", 100, 0.0, 0.25);

	TFile* f = new TFile("digi_complete.root");
	TTree* t = (TTree*)f->Get("pndsim");

	TClonesArray* clusters = new TClonesArray("PndEmcCluster");
	TClonesArray* digis = new TClonesArray("PndEmcDigi");

	t->SetBranchAddress("EmcCluster", &clusters);
	t->SetBranchAddress("EmcDigi", &digis);

	Int_t n = t->GetEntriesFast();

	PndEmcMapper::Init(1);
	for (Int_t i = 0; i < n; i++) {
		t->GetEntry(i);
		for(Int_t j = 0; j < clusters->GetEntriesFast(); j++) {
			PndEmcCluster* cluster = (PndEmcCluster*)clusters->At(j);
			PndEmcClusterEnergySums* Es = new PndEmcClusterEnergySums(*cluster, digis);
			h1->Fill(Es->E1());
			h2->Fill(Es->E9());
			h3->Fill(Es->E25());
		}
	}

	h2->SetLineColor(kRed);
	h1->SetLineColor(kBlue);
	h3->SetLineColor(kBlack);

	gStyle->SetOptTitle(0);
	gStyle->SetStatX(0.36);
	gStyle->SetStatY(0.88);
	gStyle->SetOptStat(0);

	SetStyle();
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);
	h3->SetLineWidth(2);

	h3->GetXaxis()->SetTitle("Energy(GeV)");
	h3->GetYaxis()->SetTitle("Entries");
	FormatAxis(h3->GetXaxis());
	FormatAxis(h3->GetYaxis());

	TLegend * leg = new TLegend(0.8,0.75,0.89,0.89);
	leg->AddEntry(h1, "E1", "L");
	leg->AddEntry(h2, "E9", "L");
	leg->AddEntry(h3, "E25", "L");
	leg->SetTextFont(42);

	h3->Draw();
	h2->Draw("SAME");
	h1->Draw("SAME");
	leg->Draw("SAME");
}

