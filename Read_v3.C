const int step_distance(100),step_angle(300);
void MyLoop(int CUT){
    int bin1(100),bin2(200),bin3(50);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(90),zmin(0),zmax(1.01);
    
    string file_name("doc/Shower_hit_90.txt");
    ifstream file;
    file.open(file_name, ios::in);
    
    ofstream wirtefile;
    wirtefile.open( "doc/Fit_par.txt", ios::app);
    
    int cunt(0), cunt0(0);
    
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
    //gStyle->SetOptFit(1111);
    
    //TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,2*ymax,bin3,zmin,zmax);
    TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",200,0,3,200,0,1);
    //TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",85,0,3,85,0,90);
    h2D->GetYaxis()->SetTitle("angle");
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 1);
    
    /*TF1 *f=new TF1("f","([0]*x*x*x + [1]*x*x + [2]*x + [3]) / (x*x*x + [4]*x*x + [5]*x + [6])",0,3);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(-2.337, 17.14, -42.43, 36.76 , 16.2, -45.87, 42.45);*/
    /*f->SetParLimits(0, -6.318, 1.644);
    f->SetParLimits(1, -11.67, 45.95);
    f->SetParLimits(2, -112.7, 27.83);
    f->SetParLimits(3, -22.82, 96.34);
    f->SetParLimits(4, -13.43, 45.84);
    f->SetParLimits(5, -123.1, 31.33);
    f->SetParLimits(6, -26.5, 111.4);*/
    //f->SetParLimits(29,-10,10);
    
    /*TF1 *f=new TF1("f","([0]*x*x*x + [1]*x*x + [2]*x + [3]) / (x*x*x*x*x + [4]*x*x*x*x + [5]*x*x*x + [6]*x*x + [7]*x +[8])",0,3);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(0.3556, -1.799, 3.054, -0.7111, -1.937, 2.099, -2.575, 3.595, -0.8281);*/
    
    /*TF1 *f=new TF1("f","[0] / (x*x*x*x*x + [1]*x*x*x*x + [2]*x*x*x + [3]*x*x + [4]*x +[5])",0,3);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.269, -1.524, 0.9138, 0.02971, 0.01683, 1.474);*/
    
    TF1 *f=new TF1("f","1 / ([0]*x*x*x*x*x - [1]*(x*x*x*x - 0.865*x*x*x + 0.349*x*x) + 0.241*x +1.133)",0,3);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.049,2);
    
    /*TF1 *f=new TF1("f","1 / ((0.23*[0]+0.564)*x*x*x*x*x - [0]*(x*x*x*x - 0.865*x*x*x + 0.349*x*x) + 0.241*x +1.133)",0,3);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameter(0,2);*/
    
    /*TF1 *f=new TF1("f","1 / ((1.585-0.0281*y)*x*x*x*x*x - (73.41/(y+24.672))*(x*x*x*x - 0.865*x*x*x + 0.349*x*x) + 0.241*x +1.133)",0,3);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);
    f->SetParameters(1.049,2);*/
    
    double d_distance = (xmax - xmin)/step_distance;
    double d_angle = (ymax - ymin)/step_angle;
    map<pair<int, int>,pair<int,double>> data;
    
    double dis[step_distance];
    double ang[step_angle];
    int count[step_distance][step_angle];
    double Energy[step_distance][step_angle];
    memset(count,0,sizeof(count));
    memset(count,0,sizeof(count));
    
    for (int i=0; i < step_distance; i++)
        dis[i] = (d_distance/2) + i*d_distance;
    for (int j=0; j < step_angle; j++)
        ang[j] = (d_angle/2) + j*d_angle;
    
    string str1,str2,str3;
    while (!file.eof()) {
        cunt0++;
        getline(file,str1);
        getline(file,str2);
        getline(file,str3);
        double distance= atof(str1.c_str());
        double angle= atof(str2.c_str());
        double E= atof(str3.c_str());
        int N_distance = (int)(distance/d_distance);
        int N_angle = (int)(angle/d_angle);
        if ((N_distance>=step_distance)||(N_angle>=step_angle)) continue;
        if (E<0||E>16) continue;
        count[N_distance][N_angle]++;
        Energy[N_distance][N_angle]+=E;
        cunt++;
    }
    
    Double_t rr(0.9);
    int j=CUT;
    for (int i=0; i < step_distance; i++){
        //for (int j=0; j < step_angle; j++){
            //if (j != CUT) continue;
            if (count[i][j]<=0||isnan(Energy[i][j]/count[i][j])||abs(Energy[i][j])>100) {
                //h2D->Fill(dis[i],ang[j],rr);
                h2D->Fill(dis[i],rr);
                continue;
            }
                //h2D->Fill(dis[i],ang[j],Energy[i][j]/count[i][j]);
                h2D->Fill(dis[i],Energy[i][j]/count[i][j]);
                rr = Energy[i][j]/count[i][j];
        //}
    }
    
    cout << "Passed:" << cunt << " / " << cunt0 << endl;
    c1->cd();
    //h2D->Draw("CONT4Z");
    h2D->Draw("SCAT");
    h2D->Fit(f,"R");
    //f->Draw("SAME");
    
    wirtefile << CUT << "   " << f->GetParameter(0) << "   " << f->GetParameter(1) << endl;
    
    file.close();
    wirtefile.close();
}
int Read_v3(int j = 0){
    //for (int j=0; j < step_angle; j++){
        MyLoop(j);
    //}
    return 0;
}
