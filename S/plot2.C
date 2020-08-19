void plot2(){
	TFile* f=new TFile("../../data/evtcomplete_sim.root");
	TTree* t=(TTree*)f->Get("pndsim");
	TClonesArray* hits=new TClonesArray("PndEmcHit");
	t->SetBranchAddress("EmcHit",&hits);
	
	ofstream out("out.txt",ios::app);
	Int_t n=t->GetEntries();
	for (Int_t i=0;i<n;i++){
		t->GetEntry(i);
		for(Int_t j=0;j<hits->GetEntriesFast();j++){
			PndEmcHit* hit = (PndEmcHit*)hits->At(j);
            FairMultiLinkedData flink = hit -> GetTrackEntering();
            hit->Print();
			//flink.PrintLinkInfo(out);
			out<<endl;
//			hit -> GetLinks();
//			FairLink link  = hit->GetLink(1);
		}
	}
	out.close();
//	TCanvas* c1=new TCanvas();
}
