int Read(){
    int bin1(50),bin2(50),bin3(50);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(46),zmin(0),zmax(1.01);
    string file_name("doc/Shower_hit.txt");
    //string file_name("doc/Shower_hit_90.txt");
    ifstream file;
    file.open(file_name, ios::in);
    int cunt(0);
    
    TCanvas* c1=new TCanvas("PANDA1","Hit1",tx,ty);
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
    gStyle->SetOptFit(1111);
    
    //TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,ymax,bin3,zmin,zmax);
    TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",200,0,5,200,0,1);
    h2D->GetYaxis()->SetTitle("angle");
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 1);
    
    TF1 *f=new TF1("f","[0]*exp(-[1]*pow(x,[2]))",1,2.5);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(0.5,0.5);
    //f->SetParLimits(0, -100, 100);
    //f->SetParLimits(1, -100, 100);
    //->SetParameters(2,-10,10);
    //f->SetParameters(3,-10,10);
    //f->SetParameters(4,-10,10);
    
    string str1,str2,str3;
    while (!file.eof()) {
        getline(file,str1);
        getline(file,str2);
        getline(file,str3);
        double distance= atof(str1.c_str());
        double angle= atof(str2.c_str());
        double E= atof(str3.c_str());
        if (angle>0 && angle<1)
            h2D->Fill(distance,E);
        cunt++;
    }
    cout << "Passed:" << cunt << endl;
    c1->cd();
    h2D->Draw("SCAT");
    h2D->Fit(f,"R");
    f->Draw("SAME");
    file.close();
    return 0;
}
