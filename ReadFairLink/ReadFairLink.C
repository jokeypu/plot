int ReadFairLink()
{
  FairRunAna *fRun = new FairRunAna();
  TFile* file = new TFile("../../data/evtcomplete_sim.root");
  FairFileSource* source = new FairFileSource(file,"InputFile");

  FairRootManager* ioman = FairRootManager::Instance();
  ioman->SetSource(source);
  ioman->InitSource();
//TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
  
  TClonesArray* fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
  if (!fPointArray) { cout<<"no EmcPoint"<<endl; return -1;}
  TClonesArray* fTruthArray = (TClonesArray*) ioman->GetObject("MCTrack");
  if (!fTruthArray) { cout<<"no MCTrack"<<endl; return -2;}


  int maxEvtNo = ioman->CheckMaxEventNo();
  cout<<"maxEvtNo "<<maxEvtNo<<endl;

//TH2D* h2E = new TH2D("h2E","E#gamma", 72, 0, 72, 160, 0, 160);
//h2E->GetXaxis()->SetTitle("Row");
//h2E->GetYaxis()->SetTitle("cpy&cry");
//TH2D* h2E2 = new TH2D("h2E2","E#gamma", 65, 20, 150, 180, -180, 180);
//h2E2->GetXaxis()->SetTitle("#theta");
//h2E2->GetYaxis()->SetTitle("#phi");
//
 // double ene=1.0;
//TH1D* htotE = new TH1D("htotE","E#gamma in points", 100, 0.8*ene, 1.2*ene);
//htotE->GetXaxis()->SetTitle("E_{#gamma} (in points)");
//TH1D* htotE2 = new TH1D("htotE2","E#gamma in hits", 100, 0.8*ene, 1.2*ene);
//htotE2->GetXaxis()->SetTitle("E_{#gamma} (in hits)");
//htotE2->SetLineColor(kRed);

//TH2D* h2E3 = new TH2D("h2E3","E#gamma", 65, 20, 150, 180, -180, 180);
//h2E3->GetXaxis()->SetTitle("#theta");
//h2E3->GetYaxis()->SetTitle("#phi");
    
  ofstream out("out1.txt",ios::app);
  for (int ievt=0; ievt<1; ievt++) {
     ioman->ReadEvent(ievt); // read event by event
     // read by time stamp????
     //ioman->ReadEvent(0);
     //FairEventHeader* feh = new FairEventHeader();
     //ioman->FillEventHeader(feh);
     //cout<<"Run id "<<feh->GetRunId()<<endl;

     //ioman->ReadNextEvent(1); // does not works
	 
     int npoints = fPointArray->GetEntriesFast();
     cout<<"points is "<<npoints<<endl;
     //std::map<int, std::vector<int> > TrackPoint;
	 for (int i=0; i<npoints; i++) {
       PndEmcPoint* point = (PndEmcPoint*)fPointArray->At(i);
       //point->GetLinksWithType(ioman->GetBranchId("MCTrack"));
	   //point->GetSortedMCTracks();
         
	   std::set<FairLink> links = point->GetLinks();
       out << "[ ";
       for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter;}
       out << " ]" << std::endl;
         
         
         
         
	  // int idx = links.begin()->GetIndex();
	  // if (links.size()>1) cout<<"more than 1 tracks"<<endl;
	  // if (TrackPoint.find(idx)==TrackPoint.end()) {
	 //    std::vector<int> vect;
	//	 vect.clear();
//TrackPoint[idx] = vect;
	   }
	//   TrackPoint[idx].push_back(i);

	 //cout<<i<<" "<<links.size()
	 //<<" "<<ioman->GetBranchName(links.begin()->GetType())
	 //<<" "<<links.begin()->GetIndex()
	 //<<endl;
     //vector<FairLink> links = point->GetSortedMCTracks();
	 //cout<<i<<" "<<links.size()<<endl;
	 //TVector3 pos;
	 //point->Position(pos);
	 //double eneloss = point->GetEnergyLoss();
	 //h2E3->Fill(pos.Theta()*TMath::RadToDeg(),pos.Phi()*TMath::RadToDeg(),eneloss);
     }
	 
 /*    int ntrack = fTruthArray->GetEntriesFast();
	 std::map<int, std::vector<int> > TrackDaug;
	 for (int i=0; i<ntrack; i++) {
       PndMCTrack* track  = (PndMCTrack*)fTruthArray->At(i);
	   //TVector3 pos = track->GetStartVertex();
	   if (track->GetMotherID()>=0) {
		   int idx = track->GetMotherID();
		   if (TrackDaug.find(idx)==TrackDaug.end()) {
		     std::vector<int> vect;
			 vect.clear();
			 TrackDaug[idx] = vect;
		   }
	       TrackDaug[idx].push_back(i);

	     //PndMCTrack* mothertrack  = (PndMCTrack*)fMCtrackArray->At(track->GetMotherID());
	     //TVector3 mpos = mothertrack->GetStartVertex();
	   }
	 //double eneloss = point->GetEnergyLoss();
     }
	 
	 //std::map<int, std::vector<int> >::iterator iter = TrackPoint.begin();
	 //for (;iter!=TrackPoint.end();iter++) {
	 for (int iter=0; iter<ntrack; iter++) {  
	   PndMCTrack *truth = (PndMCTrack*)fTruthArray->At(iter);
	   int motherid = truth->GetMotherID();
	   double energy = truth->Get4Momentum().E();
	   int type = truth->GetPdgCode();
	   cout<<iter<<" "<<type<<" "<<motherid<<" "<<energy*1e3<<" MeV, points: ";
	   double toteloss = 0;
	   for (int i=0; i<TrackPoint[iter].size();i++) {
	     double eloss = ((PndEmcPoint*)fPointArray->At(TrackPoint[iter].at(i)))->GetEnergyLoss();
	     toteloss += eloss;
		 cout<<" "<<TrackPoint[iter].at(i)<<":"<<eloss*1e3;
	   }
	   cout<<" tot dE="<<toteloss<<" GeV"<<endl;
	   double totedau = 0;
	   if (TrackDaug[iter].size()==0) cout<<" no daughter.";
	   else {cout<<" daughers: ";
	     for (int i=0; i<TrackDaug[iter].size();i++) {
	       PndMCTrack *truthdau = (PndMCTrack*)fTruthArray->At(TrackDaug[iter].at(i));
		   totedau += truthdau->Get4Momentum().E();
		   cout<<" "<<TrackDaug[iter].at(i);
	     }
	   }
	   cout<<" tot Edau "<<totedau<<" GeV."<<endl;
	   if (fabs(energy-toteloss-totedau)>0.00052 && TrackDaug[iter].size()==0) { 
	     cout<<"energy missing!!!"<<endl;
	     TVector3 pos = truth->GetStartVertex();
	     cout<<"position ("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<"), (r, z)=("<<pos.Perp()<<", "<<pos.Z()<<")"<<endl;
	   }
	   cout<<endl;
	 }

  }

//TCanvas* c1 = new TCanvas();
//h2E->Draw("colz");
//TCanvas* c2 = new TCanvas();
//h2E2->Draw("colz");
//TCanvas* c21 = new TCanvas();
//h2E3->Draw("colz");
//TCanvas* c3 = new TCanvas();
//htotE->Draw();
//htotE2->Draw("same");
*/
  out.close();
  return 0;
}
