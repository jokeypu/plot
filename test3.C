void savewmass() {
  gStyle->SetOptStat(1001);
  TCanvas* c1 = new TCanvas("wmass","w-mass",50,50,800,600);

  // TFile *f1 = new TFile("WW.root");
  //TFile *f2 = new TFile("TT.root");
  //TFile *f3 = new TFile("DY.root");


    
   TH1F *h1 = new TH1F("h1","h1",100,0,100);
   h1->Fill(40);
   h1->SetFillColor(6);

   TH1F *h2 = new TH1F("h2","h2",100,0,100);
   h2->Fill(50);
   h2->SetFillColor(4);

   TH1F *h3 = new TH1F("h3","h3",100,0,100);
   h3->Fill(60);
   h3->SetFillColor(3);

   THStack* hs= new THStack("hs","W-Mass stacked histograms");
    h1->Draw();
    gPad->Update();
   TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    h2->Draw();
    gPad->Update();
   TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");

    ps1->SetX1NDC(0.132832);
    ps1->SetY1NDC(0.773913);
    ps1->SetX2NDC(0.333333);
    ps1->SetY2NDC(0.874783);
    
    ps2->SetX1NDC(0.132832);
    ps2->SetY1NDC(0.624348);
    ps2->SetX2NDC(0.333333);
    ps2->SetY2NDC(0.725217);


//   c1->cd(1);
   hs->Add(h1);
  hs->Add(h2);
  hs->Add(h3);

   hs->Draw();
   ps1->Draw();
   ps2->Draw();
   gPad->Update();
  hs->GetXaxis()->SetTitle("GeV");
  hs->GetYaxis()->SetTitle("# of entries");

  TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);

  leg->AddEntry(h1,"WW","f");
  leg->AddEntry(h2,"TT","f");
  leg->AddEntry(h3,"DY","f");
  leg->Draw();
}
int test3(){
    savewmass();
    return 0;
}
