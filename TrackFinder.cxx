//PuQing
//2020.4.17
//#include "TrackFinder.h"

#include "PndMCTrack.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>

TrackFinder::TrackFinder( TClonesArray *fMCtrackArray)
:fIsSameVertex(false), mid(-1), fIsmother(2)
{
    fMCtrack = fMCtrackArray;
    TrackID.clear();
}

TrackFinder::TrackFinder( TClonesArray *fMCtrackArray, Int_t id1x, Int_t id2x)
:fIsSameVertex(false), mid(-1), id1(id1x), id2(id2x), fIsmother(2)
{
    fMCtrack = fMCtrackArray;
    TrackID.clear();
    Exec();
}

TrackFinder::~TrackFinder() {}

void TrackFinder::Print(Int_t idx){
    TVector3 svertex,smomentum;
    PndMCTrack* track = (PndMCTrack*)fMCtrack->At(idx);
    svertex = track->GetStartVertex();
    smomentum = track->GetMomentum();
    TLorentzVector s4momentum = track->Get4Momentum();
    std::cout << "StartVertex:" << "("<< svertex.x() << "," << svertex.y() << "," << svertex.z() << ")\t" << "mag:" << svertex.Mag();
    std::cout << "\tE: " << s4momentum.E();
    std::cout << "\ttheta,phi:" << "(" << smomentum.Theta() << "," << smomentum.Phi() << ")" << "\tPDG:" <<  track->GetPdgCode() << std::endl;
    
}

void TrackFinder::Print(){
    if (fIsSameVertex) { std::cout<< "id:" << id1 << " " << id2 << "\t" << "mid:" << mid << "\t";PrintIsMother(); Print(mid);}
    else{
        std::cout << "No Mother Track! " << std::endl;
    }
}

void TrackFinder::PrintPath(){
    std::cout << "\nid:" << id1 << " & " << id2 <<std::endl;
    for (Int_t i = Path1.size()-1; i > 0; i--) std::cout << Path1[i] << " -> " ;
    std::cout << Path1[0] << std::endl;
    for (Int_t i = Path2.size()-1; i > 0; i--) std::cout << Path2[i] << " -> " ;
    std::cout << Path2[0] << std::endl;
}

Int_t TrackFinder::GetPdgCode(Int_t idx){
    PndMCTrack* track = (PndMCTrack*)fMCtrack->At(idx);
    return track->GetPdgCode();
}

void TrackFinder::PrintPDG(Int_t idx){
    std::cout << "PDG(id:" << idx << ") " << GetPdgCode(idx) << "\t";
}

void TrackFinder::PrintPDG(std::vector <int> vect){
    for (int i = 0; i < vect.size();i++ ) {
        PrintPDG(vect[i]);
    }
}

void TrackFinder::IsMother(){
    fIsmother = WhoIsMother(id1, id2);
}

Int_t TrackFinder::WhoIsMother(Int_t id1x, Int_t id2x){
    Int_t IM;
    PndMCTrack* Mother1 = (PndMCTrack*)fMCtrack->At(id1x);
    Int_t Mid1 = Mother1->GetMotherID();
    PndMCTrack* Mother2 = (PndMCTrack*)fMCtrack->At(id2x);
    Int_t Mid2 = Mother2->GetMotherID();
    if ( Mid2==id1 ) IM = 1;
    else if ( Mid1==id2 ) IM = -1;
    else IM = 0;
    return IM;
}

void TrackFinder::PrintIsMother(){
    if ( fIsmother == 1 ) std::cout<< id1 << " -> " << id2 << "\t";
    else if ( fIsmother == -1 ) std::cout<< id2 << " -> " << id1 << "\t";
    else if ( fIsmother == 0 ) std::cout<< id1 << " xx " << id2 << "\t";
    else std::cout<< "-E\t";
}

void TrackFinder::AddTrackID(Int_t id1x, Int_t id2x){
    if ( id1x >= 0 ) {TrackID.push_back(id1x); id1 = id1x;}
    if ( id2x >= 0 ) {TrackID.push_back(id2x); id2 = id2x;}
    Exec();
    IsMother();
}

void TrackFinder::Exec(){
    Path1.clear();
    Path2.clear();
    Int_t mid1 = id1;
    while ( mid1 >=0 && fIsSameVertex == 0 ){
        Int_t mid2 = id2;
        Path2.clear();
        while ( mid2 >=0 ){
            if ( mid1 == mid2 ) {
                fIsSameVertex = 1;
                mid = mid1;
                Path2.push_back(mid);
                break;
            }
            Path2.push_back(mid2);
            PndMCTrack* mtrack2 = (PndMCTrack*)fMCtrack->At(mid2);
            mid2 = mtrack2->GetMotherID();
        }
        Path1.push_back(mid1);
        PndMCTrack* mtrack1 = (PndMCTrack*)fMCtrack->At(mid1);
        mid1 = mtrack1->GetMotherID();
    }
    if (!fIsSameVertex) {
        Path1.clear();
        Path2.clear();
    }
}
ClassImp(TrackFinder)
