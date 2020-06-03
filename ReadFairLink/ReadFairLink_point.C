int ReadFairLink_point()
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
  ofstream out("point_links.txt",ios::app);

  for (int ievt=0; ievt<1; ievt++) {
     ioman->ReadEvent(ievt); // read event by event
     int npoints = fPointArray->GetEntriesFast();
     cout<<"points is "<<npoints<<endl;
	 for (int i=0; i<npoints; i++) {
       PndEmcPoint* point = (PndEmcPoint*)fPointArray->At(i);

       std::set<FairLink> links = point->GetLinks();
       out << "[ ";
       for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter;}
       out << " ]" << std::endl;
	   }
     }
  out.close();
  return 0;
}
