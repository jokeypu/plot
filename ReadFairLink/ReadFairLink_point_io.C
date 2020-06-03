int ReadFairLink_point_io()
{
  FairRunAna *fRun = new FairRunAna();
  TFile* file = new TFile("../../data/evtcomplete_sim.root");
  FairFileSource* source = new FairFileSource(file,"InputFile");

  FairRootManager* ioman = FairRootManager::Instance();
  ioman->SetSource(source);
  ioman->InitSource();

  TClonesArray* fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
  if (!fPointArray) { cout<<"no EmcPoint"<<endl; return -1;}
  int maxEvtNo = ioman->CheckMaxEventNo();
  cout<<"maxEvtNo "<<maxEvtNo<<endl;
    
  //FairLinkManager* linkman = FairLinkManager::Instance();
  //std::cout << "EmcPoint Link Manager:" << linkman->IsIgnoreType(0) << std::endl;
  ofstream out("point_links.txt",ios::app);
  for (int ievt=0; ievt<1; ievt++) {
     ioman->ReadEvent(ievt); // read event by event
	 
     int npoints = fPointArray->GetEntriesFast();
     cout<<"points is "<<npoints<<endl;
       for (int i=0; i<npoints; i++) {
       PndEmcPoint* point = (PndEmcPoint*)fPointArray->At(i);
	std::cout << "Crystal: " << point->GetCrystal();
	std::cout << "\t NLink Entering: " << point->GetEntering();
	std::cout << "\t NLink Exiting: " << point->GetExiting() << std::endl;
	   }
     }
  out.close();
  return 0;
}
