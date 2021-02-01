int z_Pi_mom(){
    const double m_pi0 = 0.1349768;
    
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
    
    TH2D* h = new TH2D("Hist","h",500,1,6,500,2,15);
    //TH2D* h = new TH2D("Hist","h",500,0,6,500,2,181);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 0.5);
    h->GetXaxis()->SetTitle("E_{#pi0}   [GeV]");
    h->GetYaxis()->SetTitle("angle   [deg]");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    
    TH2D* h_E = new TH2D("Hist_E","h_E",500,0,6,500,0,6);
    h_E->SetMarkerStyle(7);
    h_E->SetMarkerColorAlpha(kAzure+3, 0.5);
    h_E->GetXaxis()->SetTitle("E_{#pi0}   [GeV]");
    h_E->GetYaxis()->SetTitle("E_{#gamma}   [GeV]");
    h_E->GetXaxis()->CenterTitle();
    h_E->GetYaxis()->CenterTitle();
    h_E->GetZaxis()->CenterTitle();
    
    for (double E_pi0 = 0; E_pi0 < 6 ; E_pi0+= 0.001){
        for (int i = 0; i < 1000 ; i++){
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
            h_E->Fill(E_pi0,E_gamma1);
            //h_E->Fill(E_pi0,E_gamma2);
            h->Fill(E_pi0,delta_angle);
        }
    }
    for (int i = 1; i < 501; i++){
        for (int j = 1; j < 501; j++){
            if (h->GetBinContent(i,j) != 0) h->SetBinContent(i, j, (int)(TMath::Log(h->GetBinContent(i,j))));
            if (h_E->GetBinContent(i,j) != 0) h_E->SetBinContent(i, j, (int)(TMath::Log(h_E->GetBinContent(i,j))));
        }
    }
    c1->cd();
    h->Draw("PCOLZ");
    
    c2->cd();
    h_E->Draw("PCOLZ");
    return 0;
}
