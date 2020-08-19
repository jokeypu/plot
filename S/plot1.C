void plot1(){
	TH1D* h1=new TH1D("h1","htemp1",100,0.,0.01);
	TH1D* h2=new TH1D("h2","htemp2",100,0.,0.01);

	TFile* f=new TFile("../../data/evtcomplete_sim.root");
	TTree* t=(TTree*)f->Get("pndsim");

	TClonesArray* hits=new TClonesArray("PndEmcHit");
	TClonesArray* points=new TClonesArray("PndEmcPoint");

	t->SetBranchAddress("EmcHit",&hits);
	t->SetBranchAddress("EmcPoint",&points);
	
	Int_t n=t->GetEntries();
	for (Int_t i=0;i<n;i++){
		t->GetEntry(i);
		for(Int_t j=0;j<hits->GetEntriesFast();j++){
			PndEmcHit* hit=(PndEmcHit*)hits->At(j);
			h1->Fill(hit->GetEnergy());
		}

		for(Int_t j=0;j<points->GetEntriesFast();j++){
			FairMCPoint* point=(FairMCPoint*)points->At(j);
			h2->Fill(point->GetEnergyLoss());
		}
	}
	TCanvas* c1=new TCanvas();
	h1->SetLineColor(kRed);
	h2->SetLineColor(kBlue);
	c1->Divide(1,2);
	c1->cd(1);
	h1->Draw();
	c1->cd(2);
	h2->Draw();
}
