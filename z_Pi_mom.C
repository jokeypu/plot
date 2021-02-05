int z_Pi_mom(){
    const double m_pi0 = 0.1349768;
    int nbins1 = 25000, nbins2 = 3000;
    Int_t BinCut =  10000, NN = 10000;
    
    TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
    TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1,"xyz");
    //gStyle->SetPalette(62);
    //gStyle->SetPalette(1);
    
    /*const Int_t Number = 3;
    Double_t Red[Number]    = { 1.00, 0.00, 0.00};
    Double_t Green[Number]  = { 0.00, 1.00, 0.00};
    Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
    Double_t Length[Number] = { 0.00, 0.50, 1.00 };
    Int_t nb=50;
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);*/
    
    
    const Int_t NColor = 200;
    Int_t colors[NColor];
    for (int i = 0; i < 49 ; i++){
        colors[i] = 51+i;
    }
    
    for (int i = 49; i < NColor ; i++){
        if (i >=49 && i<80) colors[i] = TColor::GetColor(255,  0,  0);
        if (i >=80 && i<110) colors[i] = TColor::GetColor(230,  0,  0);
        if (i >=110 && i<140) colors[i] = TColor::GetColor(205,  0,  0);
        if (i >=140 && i<170) colors[i] = TColor::GetColor(180,  0,  0);
        if (i >=170 && i<200) colors[i] = TColor::GetColor(155,  0,  0);
    }
    
    gStyle->SetPalette(NColor,colors,1);
    
    TH2D* h = new TH2D("Hist","h",nbins1,0,25,nbins1,0.5,50);
    //TH2D* h = new TH2D("Hist","h",500,0,6,500,2,181);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 0.5);
    h->GetXaxis()->SetTitle("P_{#pi^{0}}   [GeV/c]");
    h->GetYaxis()->SetTitle("angle   [deg]");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    
    TH2D* h_E = new TH2D("Hist_E","h_E",nbins2,0,30,nbins2,0,30);
    h_E->SetMarkerStyle(7);
    h_E->SetMarkerColorAlpha(kAzure+3, 0.5);
    h_E->GetXaxis()->SetTitle("P_{#pi^{0}}   [GeV/c]");
    h_E->GetYaxis()->SetTitle("E_{#gamma}   [GeV]");
    h_E->GetXaxis()->CenterTitle();
    h_E->GetYaxis()->CenterTitle();
    h_E->GetZaxis()->CenterTitle();
    
    for (double E_pi0 = 0.0005; E_pi0 < 25 ; E_pi0+= 0.001){
        for (int i = 0; i < NN ; i++){
            double rd1 = 2*(rand()/(RAND_MAX+1.))-1;
            double rd2 = 2*(rand()/(RAND_MAX+1.))-1;
            double rd3 = 2*(rand()/(RAND_MAX+1.))-1;
            TVector3 vv(rd1,rd2,rd3);
            vv.SetMag(m_pi0/2);
            double px = vv.X(), py = vv.Y(), pz = vv.Z();
            TVector3 vv_n1(-0.5*sqrt(E_pi0*E_pi0-m_pi0*m_pi0)+px*E_pi0/m_pi0,py,pz);
            TVector3 vv_n2(-0.5*sqrt(E_pi0*E_pi0-m_pi0*m_pi0)-px*E_pi0/m_pi0,-py,-pz);
            double delta_angle = TMath::RadToDeg()*(vv_n1.Angle(vv_n2));
            double E_gamma1 = E_pi0/2-px/m_pi0*sqrt(E_pi0*E_pi0-m_pi0*m_pi0);
            double E_gamma2 = E_pi0/2+px/m_pi0*sqrt(E_pi0*E_pi0-m_pi0*m_pi0);
            //cout << E_pi0 << "   " << E_gamma1 << ", " << E_gamma2 << endl;
            //h_E->Fill(E_pi0-m_pi0,E_gamma1);
            //h_E->Fill(E_pi0,E_gamma2);
            h->Fill(E_pi0-m_pi0,delta_angle);
        }
    }
    
    //Int_t BinCut =  20;
    for (int i = 1; i < nbins1+1; i++){
        for (int j = 1; j < nbins1+1; j++){
            //if (h->GetBinContent(i,j) != 0) h->SetBinContent(i, j, (int)(TMath::Log10(h->GetBinContent(i,j))));
            //if (h->GetBinContent(i,j) > 3*BinCut) h->SetBinContent(i, j, 3*BinCut);
            Double wx = h->ProfileX()->GetBinWidth(i);
            Double wy = h->ProfileY()->GetBinWidth(j)
            h->SetBinContent(i, j, h->GetBinContent(i,j)/NN/(wx*xy));
            //if (h_E->GetBinContent(i,j) > BinCut) h_E->SetBinContent(i, j, BinCut);
            //if (h_E->GetBinContent(i,j) != 0) h_E->SetBinContent(i, j, (int)(TMath::Log10(h_E->GetBinContent(i,j))));
        }
    }
    
    c1->cd();
    c1->SetGridy();
    c1->SetLogy();
    h->Draw("PCOLZ");
    
    c2->cd();
    h_E->Draw("PCOLZ");
    return 0;
}
