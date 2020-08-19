//PuQing
//2020.4.16
#ifndef TRACKFINDER_H
#define TRACKFINDER_H

#include <vector>

class TClonesArray;
class PndMCTrack;

class TrackFinder {
public:
    TrackFinder( TClonesArray *fMCtrack);
    TrackFinder( TClonesArray *fMCtrack, Int_t id1x, Int_t id2x);
    virtual ~TrackFinder();
    virtual void Print(Int_t idx);
    virtual void Print();
    void PrintPDG(Int_t idx);
    void PrintPDG(std::vector <int> vect);
    virtual void AddTrackID(Int_t idx) { TrackID.push_back(idx); };
    virtual void AddTrackID(Int_t id1x, Int_t id2x);
    bool IsSameVertex() { return fIsSameVertex; };
    void PrintIsMother();
    void PrintPath();

    Int_t GetNTracks() { return TrackID.size(); };
    Int_t WhoIsMother(Int_t id1x, Int_t id2x);
    Int_t GetPdgCode(Int_t idx);
    std::vector <int> GetPath1() { return Path1; }
    std::vector <int> GetPath2() { return Path2; }
    
private:
    std::vector <int> TrackID;
    std::vector <int> Path1;
    std::vector <int> Path2;
    std::vector <int> MTrackID;
    Int_t mid, id1, id2;
    Int_t fIsmother;
    
    PndMCTrack* track1;
    PndMCTrack* track2;
    TClonesArray *fMCtrack;
    virtual void Exec();
    void IsMother();
    
    bool fIsSameVertex;
    ClassDef(TrackFinder, 1);
};
#endif
#include "TrackFinder.cxx"
