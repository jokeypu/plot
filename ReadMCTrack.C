int ReadMCTrack()
{
	FairRunAna *fRun = new FairRunAna();
	TFile* file = new TFile("../../data/evtcomplete_sim.root");
	FairFileSource* source = new FairFileSource(file,"InputFile");
	FairRootManager* ioman = FairRootManager::Instance();
	ioman->SetSource(source);
	ioman->InitSource();

	//PndEmcMapper a;
	//PndEmcMapper::Init(1);
	//PndEmcMapper* fMapper=PndEmcMapper::Instance();
	//PndEmcTwoCoordIndex* tmpTCI = fMapper->GetTCI(108010001);
	//cout<<"x,y,idx: "<<tmpTCI->XCoord()<<","<<tmpTCI->YCoord()<<","<<tmpTCI->Index()<<endl;

	//PndEmcGeoPar* fGeoPar = (PndEmcGeoPar*) rtdb->getContainer("PndEmcGeoPar");
	//fGeoPar->InitEmcMapper();
	////par->InitEmcMapper();
	//PndEmcMapper* fMapper=PndEmcMapper::Instance();
	//PndEmcStructure* fEmcStr=PndEmcStructure::Instance();
	//const PndEmcTciXtalMap &XtalMap = fEmcStr->GetTciXtalMap(); 
	//PndEmcTwoCoordIndex *tmpTCI;
	//PndEmcXtal *tmpXtal;
	//tmpTCI = fMapper->GetTCI(108010001);
	//// tmpTCI = fMapper->GetTCI(detectorID);
	//tmpXtal =XtalMap.find(tmpTCI)->second; 
	//TVector3 frontvec = tmpXtal->frontCentre();
	//cout<<"front position "<<frontvec.Z()<<" "<<frontvec.Perp()<<endl;

	//return 0;

	TClonesArray* fMCtrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
	TClonesArray* fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
	if (!fPointArray) return -1; 
	int maxEvtNo = ioman->CheckMaxEventNo();
	cout<<"maxEvtNo "<<maxEvtNo<<endl;

	double r2d = TMath::RadToDeg();
	//TH2D* h2E1 = new TH2D("h2E2","E#gamma", 6500, 20, 150, 18000, -180, 180);
	TH2D* h2E1 = new TH2D("h2E1","position", 2500, 60, 85, 1400, 94, 108);
	h2E1->GetXaxis()->SetTitle("#theta");
	h2E1->GetYaxis()->SetTitle("#phi");
	h2E1->GetXaxis()->CenterTitle();
	h2E1->GetYaxis()->CenterTitle();
	h2E1->GetXaxis()->SetTitleSize(0.06);
	h2E1->GetYaxis()->SetTitleSize(0.06);
	h2E1->GetXaxis()->SetLabelSize(0.06);
	h2E1->GetYaxis()->SetLabelSize(0.06);

	double ene=15.0;

	TH2D* h2E3 = new TH2D("h2E3","position", 65, 20, 150, 180, -180, 180);
	h2E3->GetXaxis()->SetTitle("#theta");
	h2E3->GetYaxis()->SetTitle("#phi");

	TH2D* h2E2 = new TH2D("h2E2","postion", 100, -70, 140, 100, 0, 100);
	h2E2->GetXaxis()->SetTitle("z");
	h2E2->GetYaxis()->SetTitle("r");

	TH2D* h2E4 = new TH2D("h2E4","postion", 100, -70, 140, 100, 0, 100);
	h2E4->GetXaxis()->SetTitle("z");
	h2E4->GetYaxis()->SetTitle("r");
	h2E4->GetXaxis()->CenterTitle();
	h2E4->GetYaxis()->CenterTitle();
	h2E4->GetXaxis()->SetTitleSize(0.06);
	h2E4->GetYaxis()->SetTitleSize(0.06);
	h2E4->GetXaxis()->SetLabelSize(0.06);
	h2E4->GetYaxis()->SetLabelSize(0.06);

	gStyle->SetOptStat(0);
	TCanvas *tmpcvs = new TCanvas();
	TCanvas *tmpcvs2 = new TCanvas();
	tmpcvs->SetMargin(0.15, 0.05, 0.15, 0.05);
	tmpcvs2->SetMargin(0.15, 0.05, 0.15, 0.05);
	for (int ievt=0; ievt<1; ievt++) {
		ioman->ReadEvent(ievt); // read event by event
		// read by time stamp????
		//ioman->ReadEvent(0);
		//FairEventHeader* feh = new FairEventHeader();
		//ioman->FillEventHeader(feh);
		//cout<<"Run id "<<feh->GetRunId()<<endl;

		//ioman->ReadNextEvent(1); // does not works
		//int nhits=fHitArray->GetEntriesFast();
		//cout<<"hits is "<<nhits<<endl;
		int npoints = fPointArray->GetEntriesFast();
		cout<<"points is "<<npoints<<endl;
		std::map<int, std::vector<int> > TrackPoint;
		for (int i=0; i<npoints; i++) {
			PndEmcPoint* point = (PndEmcPoint*)fPointArray->At(i);
			//point->GetLinksWithType(ioman->GetBranchId("MCTrack"));
			//point->GetSortedMCTracks();
			std::set<FairLink> links = point->GetLinks();
			int idx = links.begin()->GetIndex();
			if (links.size()>1) cout<<"more than 1 tracks"<<endl;
			if (TrackPoint.find(idx)==TrackPoint.end()) {
				std::vector<int> vect;
				vect.clear();
				TrackPoint[idx] = vect;
			}
			TrackPoint[idx].push_back(i);
		}

		int ntrack = fMCtrackArray->GetEntriesFast();
		std::map<int, std::vector<int> > TrackDaug;
		for (int i=0; i<ntrack; i++) {
			PndMCTrack* track  = (PndMCTrack*)fMCtrackArray->At(i);
			//TVector3 pos = track->GetStartVertex();
			if (track->GetMotherID()>=0) {
				int idx = track->GetMotherID();
				if (TrackDaug.find(idx)==TrackDaug.end()) {
					std::vector<int> vect;
					vect.clear();
					TrackDaug[idx] = vect;
				}
				TrackDaug[idx].push_back(i);
			}
		}

		// draw tracks
		//int ntrack = fMCtrackArray->GetEntriesFast();
		cout<<"ntrack "<<ntrack<<endl;
		tmpcvs->Clear();
		h2E1->Draw();
		tmpcvs2->Clear();
		h2E4->Draw();
		for (int i=0; i<ntrack; i++) {
			//if (i==0) continue;
			PndMCTrack* track  = (PndMCTrack*)fMCtrackArray->At(i);
			//cout<<"id-motherid: "<<i<<"\t "<<track->GetMotherID()<<endl;
			TVector3 pos = track->GetStartVertex();
			double energy = track->Get4Momentum().E();
			int width=-1;
			int pdg = track->GetPdgCode();
			// loop TrackPoint
			int npoint = TrackPoint[i].size();
			for (int ip=0; ip<npoint; ip++) {
				PndEmcPoint* point = (PndEmcPoint*)fPointArray->At(TrackPoint[i].at(ip));
				TVector3 ppos; point->Position(ppos);
				TLine* line=0;
				TLine* line2=0;
				if (ip == 0) { // initial track
					line = new TLine(pos.Theta()*r2d, pos.Phi()*r2d, ppos.Theta()*r2d, ppos.Phi()*r2d);
					line2 = new TLine(pos.Z(), pos.Perp(), ppos.Z(), ppos.Perp());
					width = log(energy*1e3+3);
					line->SetLineWidth(width);
					line2->SetLineWidth(width);
					//if (i<10) cout<<"for init lines width "<< width<<endl;
				}
				else if (ip>0){ // track after some interactions
					PndEmcPoint* mpoint = (PndEmcPoint*)fPointArray->At(TrackPoint[i].at(ip-1));
					TVector3 mppos; mpoint->Position(mppos);
					if ((ppos-mppos).Mag()<1e-3) continue;
					line = new TLine(mppos.Theta()*r2d, mppos.Phi()*r2d, ppos.Theta()*r2d, ppos.Phi()*r2d);
					line2 = new TLine(mppos.Z(), mppos.Perp(), ppos.Z(), ppos.Perp());
					TVector3 mpmom; mpoint->Momentum(mpmom);
					width = log(mpmom.Mag()*1e3+3);
					//width = log(energy*1e3+3);
					line->SetLineWidth(width);
					line2->SetLineWidth(width);
				}
				else continue;
				if (pdg==22) {line->SetLineColor(kBlue);line2->SetLineColor(kBlue);}
				else if (pdg==11) {line->SetLineColor(kRed);line2->SetLineColor(kRed);}
				else if (pdg==-11) {line->SetLineColor(kGreen);line2->SetLineColor(kGreen);}
				else {line->SetLineColor(kYellow);line2->SetLineColor(kYellow);}
				tmpcvs->cd();
				if (i>0) line->Draw();
				tmpcvs2->cd();
				line2->Draw();
				//if (i<10) cout<<"width "<< width<<endl;
			}
			// check Track Daughter, check if they are ...
			int ndau = TrackDaug[i].size();
			if (ndau>=1) {  // photonelectron
				bool pe = true;
				PndMCTrack* trackd1  = (PndMCTrack*)fMCtrackArray->At(TrackDaug[i].at(ndau-1));
				int pdgd1 = trackd1->GetPdgCode();
				double ed1 = trackd1->Get4Momentum().E();
				double em = energy;
				if ( (pdgd1) != (11)) pe = false;
				TVector3 posd1 = trackd1->GetStartVertex();
				TVector3 posm;
				if (npoint>0) {
					PndEmcPoint* lastpoint = (PndEmcPoint*)fPointArray->At(TrackPoint[i].at(npoint-1));
					lastpoint->Position(posm);
					TVector3 mpmom; lastpoint->Momentum(mpmom);
					width = log(mpmom.Mag()*1e3+3);
					em = mpmom.Mag();
				} else {
					posm = pos;
					width = log(energy*1e3+3);
				}
				if (ed1 > em && pe && pdg==22) { 
					if (sqrt(ed1*ed1-em*em)>0.00051) // mass of electron
					{
						TLine* line = new TLine(posm.Theta()*r2d, posm.Phi()*r2d, posd1.Theta()*r2d, posd1.Phi()*r2d);
						line->SetLineWidth(width);
						line->SetLineColor(kBlue);
						tmpcvs->cd();
						if (i>0)line->Draw();
						TLine* line2 = new TLine(posm.Z(), posm.Perp(), posd1.Z(), posd1.Perp());
						line2->SetLineWidth(width);
						line2->SetLineColor(kBlue);
						tmpcvs2->cd();
						line2->Draw();
					}
				}
			}
			if (ndau>=2 && pdg==22) {  // for pair production
				PndMCTrack* trackd1  = (PndMCTrack*)fMCtrackArray->At(TrackDaug[i].at(ndau-2));
				PndMCTrack* trackd2  = (PndMCTrack*)fMCtrackArray->At(TrackDaug[i].at(ndau-1));
				int pdgd1 = trackd1->GetPdgCode();
				int pdgd2 = trackd2->GetPdgCode();
				if ( (pdgd1*pdgd2) != (-11*11)) continue;
				TVector3 posd1 = trackd1->GetStartVertex();
				TVector3 posd2 = trackd2->GetStartVertex();
				TVector3 posm;
				if (npoint>0) {
					PndEmcPoint* lastpoint = (PndEmcPoint*)fPointArray->At(TrackPoint[i].at(npoint-1));
					lastpoint->Position(posm);
					TVector3 mpmom; lastpoint->Momentum(mpmom);
					width = log(mpmom.Mag()*1e3+3);
				} else {
					posm = pos;
					width = log(energy*1e3+3);
				}
				TLine* line = new TLine(posm.Theta()*r2d, posm.Phi()*r2d, posd1.Theta()*r2d, posd1.Phi()*r2d);
				line->SetLineWidth(width);
				line->SetLineColor(kBlue);
				tmpcvs->cd();
				if (i>0) line->Draw();
				TLine* line2 = new TLine(posm.Z(), posm.Perp(), posd1.Z(), posd1.Perp());
				line2->SetLineWidth(width);
				line2->SetLineColor(kBlue);
				tmpcvs2->cd();
				line2->Draw();
			}

			//double eneloss = point->GetEnergyLoss();
			h2E3->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
			h2E2->Fill(pos.Z(),pos.Perp());
		}

	}

	//TCanvas* c1 = new TCanvas();
	//h2E1->Draw("colz");
	TCanvas* c2 = new TCanvas();
	h2E2->Draw("colz");
	TCanvas* c21 = new TCanvas();
	h2E3->Draw("colz");
	//TCanvas* c3 = new TCanvas();
	//htotE->Draw();
	//htotE2->Draw("same");
	return 0;
}

