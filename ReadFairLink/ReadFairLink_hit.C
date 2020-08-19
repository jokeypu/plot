int ReadFairLink_hit()
{
  FairRunAna *fRun = new FairRunAna();
  TFile* file = new TFile("../../data/evtcomplete_sim.root");
  FairFileSource* source = new FairFileSource(file,"InputFile");

  FairRootManager* ioman = FairRootManager::Instance();
  ioman->SetSource(source);
  ioman->InitSource();

  TClonesArray* fHitArray = (TClonesArray*) ioman->GetObject("EmcHit");
  if (!fHitArray) { cout<<"no EmcHit"<<endl; return -1;}

  int maxEvtNo = ioman->CheckMaxEventNo();
  cout<<"maxEvtNo "<<maxEvtNo<<endl;
    
  //FairLinkManager* linkman = FairLinkManager::Instance();
  //std::cout << "EmcHit Link Manager:" << linkman->IsIgnoreType(0) << std::endl;
  ofstream out("hit_links.txt",ios::app);
  for (int ievt=0; ievt<1; ievt++) {
     ioman->ReadEvent(ievt); // read event by event
	 
     int nhits = fHitArray->GetEntriesFast();
     out<<"hits is "<<nhits<<endl;

     out << "Entering:" << endl;
       for (int i=0; i<nhits; i++) {
       PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
	std::cout << "Crystal: " << hit->GetCrystal();
	std::cout << "\t NLink Entering: " << (hit->GetTrackEntering()).GetNLinks();
	std::cout << "\t NLink Exiting: " << (hit->GetTrackExiting()).GetNLinks() << std::endl;
         
       std::set<FairLink> links = (hit->GetTrackEntering()).GetLinks();
       out << "Crystal: " << hit->GetCrystal();
       out << "\t";
       for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter << " ";}
       out << std::endl;
	   }
      
       out << "\nExiting:" << endl;
         for (int i=0; i<nhits; i++) {
         PndEmcHit* hit = (PndEmcHit*)fHitArray->At(i);
         std::set<FairLink> links = (hit->GetTrackExiting()).GetLinks();
         out << "Crystal: " << hit->GetCrystal();
         out << "\t";
         for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter << " ";}
         out << std::endl;
         }
     }
  out.close();
  return 0;
}
