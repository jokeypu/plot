#include "TrackFinder.h"

void print(const std::vector<Int_t> &isAdded) {
    cout << "(";
    for (int i = 0; i < isAdded.size() - 1 ; i++) {
        cout << isAdded[i] << ", ";
    }
    cout << isAdded[isAdded.size()] << ")" << endl;
}

void hist(TClonesArray* fDigiArray, const std::vector<Int_t> &isAdded) {
    int bin1(28),bin2(28);
    TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",bin1,0.75,1.75,bin2,0,1);
    TH2D* histxy1=new TH2D("hvx0vy01","vx vs vy1",bin1,0.75,1.75,bin2,0,1);
    TH2D* histxy2=new TH2D("hvx0vy02","vx vs vy2",bin1,0.75,1.75,bin2,0,1);
    TH2D* histxy3=new TH2D("hvx0vy03","vx vs vy3",bin1,0.75,1.75,bin2,0,1);
    TH2D* histxy4=new TH2D("hvx0vy04","vx vs vy4",bin1,0.75,1.75,bin2,0,1);
    TH2D* histxy5=new TH2D("hvx0vy05","vx vs vy5",bin1,0.75,1.75,bin2,0,1);
    for (int i = 0; i < isAdded.size(); i++) {
            PndEmcDigi* digi = (PndEmcDigi* ) fDigiArray->At(i);
            double theta = digi->GetTheta();
            double phi = digi->GetPhi();
            double EE = 8000*sqrt(sqrt((digi->GetEnergy())));
            //cout << EE << endl;
            for (int k = 0; k < EE; k++) {
                switch (isAdded[i]) {
                    case 0:
                        histxy->Fill(theta,phi);
                        break;
                    case 1:
                        histxy1->Fill(theta,phi);
                        break;
                    case 2:
                        histxy2->Fill(theta,phi);
                        break;
                    case 3:
                        histxy3->Fill(theta,phi);
                        break;
                    case 4:
                        histxy4->Fill(theta,phi);
                        break;
                    case 5:
                        histxy5->Fill(theta,phi);
                        break;
                    case -1:
                        break;
                    default:
                        cout << "-E-" << endl;
                        break;
                }
            }
    }
    TCanvas* c1=new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    histxy->GetXaxis()->SetTitle("#theta");
    histxy->GetYaxis()->SetTitle("#phi");
    histxy->SetFillColor(46);
    histxy->Draw("BOX");
    histxy1->SetFillColor(13);
    histxy1->Draw("SAMEBOX");
    histxy2->SetFillColor(7);
    histxy2->Draw("SAMEBOX");
    histxy3->SetFillColor(9);
    histxy3->Draw("SAMEBOX");
    histxy4->SetFillColor(36);
    histxy4->Draw("SAMEBOX");
    histxy5->SetFillColor(6);
    histxy5->Draw("SAMEBOX");
}

int ReadMCTrack_test6()
{
    FairRunAna *fRun = new FairRunAna();
    PndEmcMapper *emcMap = PndEmcMapper::Instance();
    TFile* file = new TFile("mydata/evtcomplete_digi.root");
    FairFileSource* source = new FairFileSource(file,"InputFile");
    FairRootManager* ioman = FairRootManager::Instance();
    ioman->SetSource(source);
    ioman->InitSource();
    
    TFile* f = new TFile("mydata/evtcomplete_sim.root");
    TTree* t = (TTree*)f->Get("pndsim");
    TClonesArray* fMCtrackArray = new TClonesArray("PndMCTrack");
    t->SetBranchAddress("MCTrack",&fMCtrackArray);
    if (!fMCtrackArray) return -1;
    
    //TH1D* h1=new TH1D("h1","htemp1",200,0.,1.2);
    //TH2D* histxy=new TH2D("hvx0vy0","vx vs vy",300,-1.6,1.6,300,50,80);
    
    TClonesArray* fClusterArray = (TClonesArray*) ioman->GetObject("EmcCluster");
    if (!fClusterArray) return -1;
    TClonesArray* fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigi");
    if (!fDigiArray) return -1;
    
    ioman->ReadEvent(0);
    Int_t clusterNr = 0;
    Double_t dt = 0;
    Double_t deltaT = 0;
    Int_t nDigisPassed = fDigiArray->GetEntriesFast();
    //std::vector<Int_t> neighbours{3, 1, 2, 3, 2, 3, 4, 0, 3, 4, 8, 9, 3, 8, 9, 10, 4, 6, 11, 12, 13, 3, 7, 12, 13, 1, 13, 3, 9, 14, 15, 2, 10, 15, 1, 15, 1, 12, 1, 13, 0, 1, 16, 0}; // make arrays dynamic, as we don't know its size in advance
    std::vector<Int_t> neighbours{1, 2, 1, 4, 3, 5, 6, 7, 3, 7, 8, 9, 1, 9, 4, 6, 11, 12, 13, 4, 7, 12, 13, 14, 4, 8, 13, 14, 15, 4, 9, 14, 15, 16, 3, 15, 16, 17, 1, 18, 1, 12, 2, 13, 21, 3, 14, 21, 22, 3, 15, 22, 23, 3, 16, 22, 23, 2, 17, 23, 0, 1, 24, 2, 20, 24, 0, 0, 1, 23, 0, 0};
    std::vector<Int_t> DigiPassed;
    std::vector<Int_t> XPadPassed;
    std::vector<Int_t> YPadPassed;
    std::vector<Double_t> TPassed;
    std::vector<Int_t> isAdded;
    std::vector<Int_t> similarities;

    for (Int_t a = 0; a < nDigisPassed; a++)
    isAdded.push_back(-1); // initialise all entries of isAdded to -1 (needed by algorithm)
    
    for (Int_t iDigi = 0; iDigi < nDigisPassed; ++iDigi)
    {
        DigiPassed.push_back(iDigi);
    }
/*
 for (Int_t d = 0; d < nDigisPassed - 1; d++)
    { // check all pairs of digis, so all digis up to the second-to-last one
        Int_t nNeighbours = 0;
        neighbours.push_back(0); // placeholder for nr of neighbours
        for (Int_t e = d + 1; e < nDigisPassed; e++)
        {
            if (FairRunAna::Instance()->IsTimeStamp())
                dt = TMath::Abs(TPassed[e] - TPassed[d]); // also take into account that non-consecutive digi pairs may differ in time more than dtau ns

            if (dt <= deltaT && static_cast<PndEmcDigi *>(fDigiArray->At(DigiPassed[e]))->isNeighbour(static_cast<PndEmcDigi *>(fDigiArray->At(DigiPassed[d]))))
            {                             // BUILT-IN PandaRoot VERSION. Check which digis are neighbours
                neighbours.push_back(e); // construct array containing neighbouring digis
                nNeighbours++;             // keep track of nr of neighbours
            }
        }
        neighbours[neighbours.size() - (nNeighbours + 1)] = nNeighbours; // write nr of neighbours to the appropiate entry in neighbours[]
    }
*/
    // Primary Clustering
    Int_t k = 0;
    Int_t l = 0;
    Int_t nClusters = 0; // nr of preclusters
    Int_t simLength = 0;
    while (l < Int_t(neighbours.size()))
    {
        if (isAdded[k] < 0)
        {                            // if it hasn't been added yet
            isAdded[k] = nClusters; // put internal cluster nr in isAdded array
            for (Int_t j = 1; j < neighbours[l] + 1; j++)
            {                                  // for the first neighbouring hit upto #neighbours for this digi
                Int_t m = neighbours[l + j]; // m = hit index
                if (isAdded[m] < 0)
                    isAdded[m] = nClusters; // if hit hasn't been added, put internal cluster nr here
                else if (isAdded[m] != nClusters)
                {
                    if (nClusters > isAdded[m])
                    {
                        similarities.push_back(nClusters);  // if it has, put current cluster nr in this array
                        similarities.push_back(isAdded[m]); // together with its cluster nr (could be different, which is why is being stored here for later comparison)
                    }
                    else
                    {
                        similarities.push_back(isAdded[m]); // different order, so we always have the largest cluster nr first
                        similarities.push_back(nClusters);
                    }
                    simLength += 2; // increment simLength by 2, as we just wrote 2 elements to similarities
                }
            }
            nClusters++;
        }
        else
        { // if isAdded[k] isn't -1, it means the digi has been added already and isAdded[k] contains its cluster nr
            for (Int_t j = 1; j < neighbours[l] + 1; j++)
            {
                Int_t m = neighbours[l + j];
                if (isAdded[m] < 0)
                    isAdded[m] = isAdded[k]; // same as before, only use cluster nr found in isAdded[k]
                else if (isAdded[m] != isAdded[k])
                {
                    if (isAdded[k] > isAdded[m])
                    {
                        similarities.push_back(isAdded[k]); // similarities[] keeps tracks of which preclusters have to be merged
                        similarities.push_back(isAdded[m]);
                    }
                    else
                    {
                        similarities.push_back(isAdded[m]); // different order, so we always have the largest cluster nr first
                        similarities.push_back(isAdded[k]);
                    }
                    simLength += 2;
                }
            }
        }
        k++;
        l += neighbours[l] + 1;
    }
    

    if (nDigisPassed != 0)
        if (isAdded[nDigisPassed - 1] < 0)
            isAdded[nDigisPassed - 1] = nClusters++;
    
    print(isAdded);
    print(similarities);
    hist(fDigiArray,isAdded);

    // Secondary clustering
    for (Int_t i = 0; i < simLength; i += 2)
    { // use info from similarities[] to merge clusters
        if (similarities[i] != similarities[i + 1])
        {
            for (Int_t m = 0; m < nDigisPassed; m++)
            {
                if (isAdded[m] == similarities[i])
                    isAdded[m] = similarities[i + 1];
            }
            for (Int_t j = i + 2; j < simLength; j++)
            {
                if (similarities[j] == similarities[i])
                {
                    similarities[j] = similarities[i + 1];
                }
                if (similarities[j+1] == similarities[i])
                {
                    similarities[j+1] = similarities[i + 1];
                }
            }
        }
    }
    print(isAdded);
    //hist(fDigiArray,isAdded);
    
    for (Int_t i = 0; i < nClusters; i++)
    {
        Int_t n = 0;
        for (Int_t j = 0; j < nDigisPassed; j++)
        {
            if (isAdded[j] == i)
            {
                n++;
                isAdded[j] = clusterNr; // info in isAdded: digi "index" belongs to cluster "isAdded[index]"
            }
        }
        if (n > 0)
            clusterNr++; // cluster nr, don't reset it while searching for clusters in a timebunch
    }
 
    print(isAdded);
    //hist(fDigiArray,isAdded);
    
    return 0;
}

