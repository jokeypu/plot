int ReadFairLink_digi()
{
  FairRunAna *fRun = new FairRunAna();
  TFile* file = new TFile("../../data/evtcomplete_digi.root");
  FairFileSource* source = new FairFileSource(file,"InputFile");

  FairRootManager* ioman = FairRootManager::Instance();
  ioman->SetSource(source);
  ioman->InitSource();

  TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
  if (!fDigiArray) { cout<<"no EmcDigi "<<endl; return -1;}

  int maxEvtNo = ioman->CheckMaxEventNo();
  cout<<"maxEvtNo "<<maxEvtNo<<endl;
    
  //FairLinkManager* linkman = FairLinkManager::Instance();
  //std::cout << "EmcDigi Link Manager:" << linkman->IsIgnoreType(0) << std::endl;
  ofstream out("digi_links.txt",ios::app);
  for (int ievt=0; ievt<1; ievt++) {
     ioman->ReadEvent(ievt); // read event by event
	 
     int ndigis = fDigiArray->GetEntriesFast();
     out << "digis is " << ndigis <<endl;
      out << "Entering:" << endl;
     // out << "\nExiting:" << endl;
       for (int i=0; i < ndigis; i++) {
       PndEmcDigi* digi = ( PndEmcDigi* ) fDigiArray->At(i);
//	std::cout << "Crystal: " << digi->GetCrystal();
	std::cout << "\t NLink Entering: " << (digi->GetTrackEntering()).GetNLinks();
	std::cout << "\t NLink Exiting: " << (digi->GetTrackExiting()).GetNLinks() << std::endl;
         
       std::set<FairLink> links = (digi->GetTrackEntering()).GetLinks();
      // std::set<FairLink> links = (digi->GetTrackExiting()).GetLinks();
       out << "\t";
       for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter << " ";}
       out << std::endl;
	   }
     }
  out.close();
  return 0;
}
