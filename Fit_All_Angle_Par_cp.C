int Fit_All_Angle_Par_cp(){
    std::string in_name = "doc/AllPar_cp.txt";
    std::ifstream par_file;
    par_file.open(in_name,std::ios::in);
    double ty0 = 1200, tx0 = 500;
    double ty1 = 833, tx1 = 500;

    TCanvas* c1=new TCanvas("PANDA0","All Angle p1",tx0,ty0);
    TCanvas* c2=new TCanvas("PANDA1","All Angle p2",tx1,ty1);
    TCanvas* c3=new TCanvas("PANDA3","All Angle p3",tx1,ty1);
    TCanvas* c4=new TCanvas("PANDA2","All Angle p4",tx0,ty0);
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
    gStyle->SetTitleOffset(1.0,"xyz");
    //gStyle->SetOptFit(1111);
    
    TGraph *g11 = new TGraph();
    g11->SetMarkerStyle(22);
    g11->SetMarkerColorAlpha(kGreen+1, 1);
    g11->GetXaxis()->SetTitle("Angle(deg)");
    g11->GetYaxis()->SetTitle("p11");
    g11->GetXaxis()->CenterTitle();
    g11->GetYaxis()->CenterTitle();
    
    TGraph *g12 = new TGraph();
    g12->SetMarkerStyle(22);
    g12->SetMarkerColorAlpha(kGreen+1, 1);
    g12->GetXaxis()->SetTitle("Angle(deg)");
    g12->GetYaxis()->SetTitle("p12");
    g12->GetXaxis()->CenterTitle();
    g12->GetYaxis()->CenterTitle();
    
    TGraph *g13 = new TGraph();
    g13->SetMarkerStyle(22);
    g13->SetMarkerColorAlpha(kRed, 1);
    g13->GetXaxis()->SetTitle("Angle(deg)");
    g13->GetYaxis()->SetTitle("p13");
    g13->GetXaxis()->CenterTitle();
    g13->GetYaxis()->CenterTitle();
    
    TGraph *g21 = new TGraph();
    g21->SetMarkerStyle(22);
    g21->SetMarkerColorAlpha(kGreen+1, 1);
    g21->GetXaxis()->SetTitle("Angle(deg)");
    g21->GetYaxis()->SetTitle("p21");
    g21->GetXaxis()->CenterTitle();
    g21->GetYaxis()->CenterTitle();
    
    TGraph *g22 = new TGraph();
    g22->SetMarkerStyle(22);
    g22->SetMarkerColorAlpha(kRed, 1);
    g22->GetXaxis()->SetTitle("Angle(deg)");
    g22->GetYaxis()->SetTitle("p22");
    g22->GetXaxis()->CenterTitle();
    g22->GetYaxis()->CenterTitle();
    
    TGraph *g31 = new TGraph();
    g31->SetMarkerStyle(22);
    g31->SetMarkerColorAlpha(kGreen+1, 1);
    g31->GetXaxis()->SetTitle("Angle(deg)");
    g31->GetYaxis()->SetTitle("p31");
    g31->GetXaxis()->CenterTitle();
    g31->GetYaxis()->CenterTitle();
    
    TGraph *g32 = new TGraph();
    g32->SetMarkerStyle(22);
    g32->SetMarkerColorAlpha(kRed, 1);
    g32->GetXaxis()->SetTitle("Angle(deg)");
    g32->GetYaxis()->SetTitle("p32");
    g32->GetXaxis()->CenterTitle();
    g32->GetYaxis()->CenterTitle();
    
    TGraph *g41 = new TGraph();
    g41->SetMarkerStyle(22);
    g41->SetMarkerColorAlpha(kGreen+1, 1);
    g41->GetXaxis()->SetTitle("Angle(deg)");
    g41->GetYaxis()->SetTitle("p41");
    g41->GetXaxis()->CenterTitle();
    g41->GetYaxis()->CenterTitle();
    
    TGraph *g42 = new TGraph();
    g42->SetMarkerStyle(22);
    g42->SetMarkerColorAlpha(kGreen+1, 1);
    g42->GetXaxis()->SetTitle("Angle(deg)");
    g42->GetYaxis()->SetTitle("p42");
    g42->GetXaxis()->CenterTitle();
    g42->GetYaxis()->CenterTitle();
    
    TGraph *g43 = new TGraph();
    g43->SetMarkerStyle(22);
    g43->SetMarkerColorAlpha(kRed, 1);
    g43->GetXaxis()->SetTitle("Angle(deg)");
    g43->GetYaxis()->SetTitle("p43");
    g43->GetXaxis()->CenterTitle();
    g43->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    while (std::getline(par_file, str)) {
        std::stringstream strStream(str);
        int NO_Angle(-1);
        float p11, p12, p13, p21, p22, p31, p32, p41, p42, p43;
        strStream >> NO_Angle >> p11 >> p12 >> p13 >> p21 >> p22 >> p31 >> p32 >> p41 >> p42 >> p43;
        
        Double_t t_min, t_max;
        if (NO_Angle == 1) {t_min = 23.8514; t_max = 24.6978;}
        if (NO_Angle == 2) {t_min = 26.4557; t_max = 27.3781;}
        if (NO_Angle == 3) {t_min = 29.4579; t_max = 30.4916;}
        if (NO_Angle == 4) {t_min = 32.6536; t_max = 33.7759;}
        if (NO_Angle == 5) {t_min = 36.1172; t_max = 37.3507;}
        if (NO_Angle == 6) {t_min = 39.9051; t_max = 41.2390;}
        if (NO_Angle == 7) {t_min = 44.2385; t_max = 45.7355;}
        if (NO_Angle == 8) {t_min = 48.8451; t_max = 50.4459;}
        if (NO_Angle == 9) {t_min = 53.7548; t_max = 55.4790;}
        if (NO_Angle == 10) {t_min = 59.0059; t_max = 60.8229;}
        if (NO_Angle == 11) {t_min = 64.7855; t_max = 66.7591;}
        if (NO_Angle == 12) {t_min = 70.8088; t_max = 72.8652;}
        if (NO_Angle == 13) {t_min = 77.0506; t_max = 79.1942;}
        if (NO_Angle == 14) {t_min = 83.4997; t_max = 85.6749;}
        if (NO_Angle == 15) {t_min = 90.2068; t_max = 92.4062;}
        if (NO_Angle == 16) {t_min = 96.8200; t_max = 99.0099;}
        if (NO_Angle == 17) {t_min = 103.361; t_max = 105.534;}
        if (NO_Angle == 18) {t_min = 109.793; t_max = 111.893;}
        if (NO_Angle == 19) {t_min = 116.067; t_max = 118.019;}
        if (NO_Angle == 20) {t_min = 121.838; t_max = 123.686;}
        if (NO_Angle == 21) {t_min = 127.273; t_max = 129.033;}
        if (NO_Angle == 22) {t_min = 132.400; t_max = 134.031;}
        if (NO_Angle == 23) {t_min = 137.230; t_max = 138.679;}
        Double_t angle = (t_max+t_min)/20.0;
        if (NO_Angle == -1 || 10*angle<18 || 10*angle > 145) return 1;
        
        g11->SetPoint(N,10*angle,p11);
        g12->SetPoint(N,10*angle,p12);
        g13->SetPoint(N,10*angle,p13);
        
        g21->SetPoint(N,10*angle,p21);
        g22->SetPoint(N,10*angle,p22);
        
        g31->SetPoint(N,10*angle,p31);
        g32->SetPoint(N,10*angle,p32);
        
        g41->SetPoint(N,10*angle,p41);
        g42->SetPoint(N,10*angle,p42);
        g43->SetPoint(N,10*angle,p43);
        
        N++;
    }
    par_file.close();

    TF1* f11=new TF1("f11","[0]",20,140);
    f11->SetLineColor(kRed-7);
    
    TF1* f12=new TF1("f12","[0]",20,140);
    f12->SetLineColor(kRed-7);
    
    TF1* f13=new TF1("f13","[0]*(x-[1])*(x-[1])+[2]",20,140);
    f13->SetLineColor(kRed-7);
    f13->SetParameters(0.04,80,0.004);
    
    TF1* f21=new TF1("f21","[0]",20,140);
    f21->SetLineColor(kRed-7);
    
    TF1* f22=new TF1("f22","[0]",20,140);
    f22->SetLineColor(kRed-7);
    
    TF1* f31=new TF1("f31","[0]",20,140);
    f31->SetLineColor(kRed-7);
    
    TF1* f32=new TF1("f32","[0]*(x-[1])*(x-[1])+[2]",20,140);
    f32->SetLineColor(kRed-7);
    f32->SetParameters(-0.04,80,0.004);
    
    TF1* f41=new TF1("f41","[0]",20,140);
    f41->SetLineColor(kRed-7);
    
    TF1* f42=new TF1("f42","[0]",20,140);
    f42->SetLineColor(kRed-7);
    
    TF1* f43=new TF1("f43","[0]*(x-[1])*(x-[1])+[2]",20,140);
    f43->SetLineColor(kRed-7);
    f43->SetParameters(0.04,80,0.004);
    
    g11->Fit(f11,"R");
    g12->Fit(f12,"R");
    g13->Fit(f13,"R");
    g21->Fit(f21,"R");
    g22->Fit(f22,"R");
    g31->Fit(f31,"R");
    g32->Fit(f32,"R");
    g41->Fit(f41,"R");
    g42->Fit(f42,"R");
    g43->Fit(f43,"R");
    
    g11->GetYaxis()->SetRangeUser(-6.5,2.5);
    g12->GetYaxis()->SetRangeUser(-4,14);
    g13->GetYaxis()->SetRangeUser(3,3.35);
    g21->GetYaxis()->SetRangeUser(-0.1,0.1);
    g22->GetYaxis()->SetRangeUser(0.6,1.3);
    g31->GetYaxis()->SetRangeUser(-0.01,0.01);
    g32->GetYaxis()->SetRangeUser(0.6,0.75);
    g41->GetYaxis()->SetRangeUser(-6,0);
    g42->GetYaxis()->SetRangeUser(0,5);
    g43->GetYaxis()->SetRangeUser(4.2,6);
    
    //c1->Divide(3, 1);
    c1->Divide(1, 3);
    c1->cd(1);
    g11->Draw("AP.");
    c1->cd(2);
    g12->Draw("AP.");
    c1->cd(3);
    g13->Draw("AP.");
    
    //c2->Divide(2, 1);
    c2->Divide(1, 2);
    c2->cd(1);
    g21->Draw("AP.");
    c2->cd(2);
    g22->Draw("AP.");
    
    //c3->Divide(2, 1);
    c3->Divide(1, 2);
    c3->cd(1);
    g31->Draw("AP.");
    c3->cd(2);
    g32->Draw("AP.");
    
    //c4->Divide(3, 1);
    c4->Divide(1, 3);
    c4->cd(1);
    g41->Draw("AP.");
    c4->cd(2);
    g42->Draw("AP.");
    c4->cd(3);
    g43->Draw("AP.");
    
    /*cout << "Double_t mp0[3] = {" << f01->GetParameter(0) << ", " << f02->GetParameter(0) << ", " << f03->GetParameter(0) << "};" << endl;
    cout << "Double_t mp1[3] = {" << f11->GetParameter(0) << ", " << f12->GetParameter(0) << ", " << f13->GetParameter(0) << "*pow(ShowerAngle-" << f13->GetParameter(1) << ",2)+" << f13->GetParameter(2) << "};" << endl;
    cout << "Double_t mp2[3] = {" << f21->GetParameter(0) << ", " << f22->GetParameter(0) << ", " << f23->GetParameter(0) << "*pow(ShowerAngle-" << f23->GetParameter(1) << ",2)+" << f23->GetParameter(2) << "};" << endl;
    cout << "Double_t mp3[2] = {" << f31->GetParameter(0) << ", " << f32->GetParameter(0) << "};" << endl;
    cout << "Double_t mp4[3] = {" << f41->GetParameter(0) << ", " << f42->GetParameter(0) << ", " << f43->GetParameter(0) << "*pow(ShowerAngle-" << f43->GetParameter(1) << ",2)+" << f43->GetParameter(2) << "};" << endl;*/
    
    //cout << "Double_t mp0[3] = {" << f01->GetParameter(0) << "*pow(ShowerAngle-" << f01->GetParameter(1) << ",2)+" << f01->GetParameter(2) << ", " << f02->GetParameter(0) << ", "
    //<< f03->GetParameter(0) << "*pow(ShowerAngle-" << f03->GetParameter(1) << ",2)+" << f03->GetParameter(2) << "};" << endl;
    //cout << "Double_t mp1[2] = {" << f11->GetParameter(0) << ", " << f12->GetParameter(0) << "};" << endl;
    //cout << "Double_t mp2[3] = {" << f21->GetParameter(0) << ", " << f22->GetParameter(0) << ", "
    //<< f23->GetParameter(0) << "*pow(ShowerAngle-" << f23->GetParameter(1) << ",2)+" << f23->GetParameter(2) << "+" << f23->GetParameter(3) << "*TMath::Gaus(ShowerAngle," << f23->GetParameter(4) << "," << f23->GetParameter(5) << ")" << "};" << endl;
    //cout << "Double_t mp3[2] = {" << f31->GetParameter(0) << ", " << f32->GetParameter(0) << "};" << endl;
    //cout << "Double_t mp4[2] = {" << f41->GetParameter(0) << ", " << f42->GetParameter(0) << "};" << endl;
    
    /*c1->Print("doc/All_Angle_FitPar_p1.png");
    c2->Print("doc/All_Angle_FitPar_p2.png");
    c3->Print("doc/All_Angle_FitPar_p3.png");
    c4->Print("doc/All_Angle_FitPar_p4.png");*/
    return 0;
}
