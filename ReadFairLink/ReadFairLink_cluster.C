int ReadFairLink_cluster()
{
  FairRunAna *fRun = new FairRunAna();
  TFile* file = new TFile("../../data/evtcomplete_digi.root");
  FairFileSource* source = new FairFileSource(file,"InputFile");

  FairRootManager* ioman = FairRootManager::Instance();
  ioman->SetSource(source);
  ioman->InitSource();

  TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
  if (!fClusterArray) { cout<<"no EmcCluster"<<endl; return -1;}

  int maxEvtNo = ioman->CheckMaxEventNo();
  cout<<"maxEvtNo "<<maxEvtNo<<endl;
    
  //FairLinkManager* linkman = FairLinkManager::Instance();
  //std::cout << "EmcCluster Link Manager:" << linkman->IsIgnoreType(0) << std::endl;
  ofstream out("cluster_links.txt",ios::app);
  for (int ievt=0; ievt<1; ievt++) {
     ioman->ReadEvent(ievt); // read event by event
	 
     int nclusters = fClusterArray->GetEntriesFast();
     out << "clusters is " << nclusters <<endl;
      out << "Entering:" << endl;
       for (int i=0; i < nclusters; i++) {
       PndEmcCluster* cluster = (PndEmcCluster*) fClusterArray->At(i);
//	std::cout << "Crystal: " << cluster->GetCrystal();
	std::cout << "\t NLink Entering: " << (cluster->GetTrackEntering()).GetNLinks();
	std::cout << "\t NLink Exiting: " << (cluster->GetTrackExiting()).GetNLinks() << std::endl;
         
       std::set<FairLink> links = (cluster->GetTrackEntering()).GetLinks();
       //std::set<FairLink> links = (cluster->GetTrackExiting()).GetLinks();
       out << "\t";
       for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter << " ";}
       out << std::endl;
	   }
      
            out << "\nExiting:" << endl;
             for (int i=0; i < nclusters; i++) {
             PndEmcCluster* cluster = (PndEmcCluster*) fClusterArray->At(i);
             std::set<FairLink> links = (cluster->GetTrackExiting()).GetLinks();
             out << "\t";
             for (std::set<FairLink>::iterator linkIter = links.begin(); linkIter != links.end(); linkIter++) { out << *linkIter << " ";}
             out << std::endl;
             }
     }
  out.close();
  return 0;
}
