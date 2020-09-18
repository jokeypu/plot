Double_t func_x0(Double_t d, Double_t par1_0, Double_t par1_1, Double_t par1_2, Double_t par1_3, Double_t par1_4, Double_t par2_0, Double_t par2_1, Double_t par2_2, Double_t par2_3, Double_t par2_4){
    if (d < 1.7) return par1_0*TMath::Vavilov(d - par1_1, par1_2, par1_3) + par1_4;
    else return par2_0*TMath::Landau(d - par2_1, par2_2, par2_3) + par2_4;
}
Double_t func_a(Double_t d, Double_t par3_0, Double_t par3_1, Double_t par3_2, Double_t par3_3, Double_t par3_4, Double_t par4_0, Double_t par4_1, Double_t par4_2, Double_t par4_3, Double_t par4_4){
    if (d < 1.39) return par3_0*TMath::Poisson(par3_1*(d-par3_2), par3_3)+par3_4;
    else return par4_0*TMath::Poisson(par4_1*(d-par4_2), par4_3)+par4_4;
}
Double_t func_h(Double_t d, Double_t par5_0, Double_t par5_1, Double_t par5_2, Double_t par6_0, Double_t par6_1, Double_t par6_2, Double_t par6_3, Double_t par7_0, Double_t par7_1, Double_t par7_2){
    if (d < 1.4 ) return par5_0*pow(d,par5_1)+par5_2;
    else if ( d < 3.5) return par6_0*TMath::Landau(d-par6_1,par6_2,par6_3);
    else return par7_0*exp(par7_1*d+par7_2);
}

Double_t m(Double_t d, Double_t angle, Double_t par1_0, Double_t par1_1, Double_t par1_2, Double_t par1_3, Double_t par1_4, Double_t par2_0, Double_t par2_1, Double_t par2_2, Double_t par2_3, Double_t par2_4, Double_t par3_0, Double_t par3_1, Double_t par3_2, Double_t par3_3, Double_t par3_4, Double_t par4_0, Double_t par4_1, Double_t par4_2, Double_t par4_3, Double_t par4_4, Double_t par5_0, Double_t par5_1, Double_t par5_2, Double_t par6_0, Double_t par6_1, Double_t par6_2, Double_t par6_3, Double_t par7_0, Double_t par7_1, Double_t par7_2){
    Double_t h = func_h(d,  par5_0,  par5_1,  par5_2,  par6_0,  par6_1,  par6_2,  par6_3,  par7_0,  par7_1,  par7_2);
    if ((d < 0.8) || (d > 2.8)) return h;
    else{
        Double_t a = func_a(d,  par3_0,  par3_1,  par3_2,  par3_3,  par3_4,  par4_0,  par4_1,  par4_2,  par4_3,  par4_4);
        Double_t x0 = func_x0(d,  par1_0,  par1_1,  par1_2,  par1_3,  par1_4,  par2_0,  par2_1,  par2_2,  par2_3,  par2_4);
        a *= a;
        angle = abs(fmod(angle,45.0) - 45*(((int)(angle/45.0))%2));
        if (angle<x0) return a*angle*angle+h;
        else return a*x0*x0+(a/(45/x0-1))*(x0-45)*(x0-45)-(a/(45/x0-1))*(angle-45)*(angle-45)+h;
    }
}

int Read_v2(){
    int bin1(100),bin2(200),bin3(50);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(90),zmin(0),zmax(1.01);
    const int step_distance(85),step_angle(85);
    
    string file_name("doc/Shower_hit_90.txt");
    ifstream file;
    file.open(file_name, ios::in);
    
    ofstream x, y, z, xyz;
    x.open("doc/X.txt",ios::out);
    y.open("doc/Y.txt",ios::out);
    z.open("doc/Z.txt",ios::out);
    xyz.open("doc/XYZ.txt",ios::out);
    
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
    gStyle->SetOptFit(1111);
    
    //TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,2*ymax,bin3,zmin,zmax);
    //TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",200,0,5,200,0,1);
    TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",85,0,3,85,0,90);
    h2D->GetYaxis()->SetTitle("angle");
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 1);
    
    //TF2 *f=new TF2("f","m(x, y, [0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14], [15], [16], [17], [18], [19], [20], [21], [22], [23], [24], [25], [26], [27], [28], [29])",0,3);
    TF2 *f=new TF2("f","m(x, y, 8.62369, 5.32678, 12.0, 0.696003, 22.8543, 83.4757, -0.171518, 1.88588, 0.199225, 23.0441, 0.0227033, 4.97094, 0.772748, 4.0, 0.0018315, 0.00790974, 1.56139, 1.39, 0.428328, 0.00115601, -0.0849149, 2.58396, 0.419646, 10.0, 0.0124098, 0.0780432, 0.182612, 0.144939, -0.435278, 0.0642399)",0,3);
    //f->SetLineWidth(2);
    //f->SetLineColor(kRed);
    //f->SetParameters(8.62369, 5.32678, 12.0, 0.696003, 22.8543, 83.4757, -0.171518, 1.88588, 0.199225, 23.0441, 0.0227033, 4.97094, 0.772748, 4.0, 0.0018315, 0.00790974, 1.56139, 1.39, 0.428328, 0.00115601, -0.0849149, 2.58396, 0.419646, 10.0, 0.0124098, 0.0780432, 0.182612, 0.144939, -0.435278, 0.0642399);
    //f->SetParameter(0,8.62369);
    /*f->SetParameters(8.62369, 5.32678, 12.0, 0.696003, 22.8543, 83.4757, -0.171518, 1.88588, 0.199225, 23.0441, 0.0227033);
    f->SetParameter(11,4.97094);
    f->SetParameter(12,0.772748);
    f->SetParameter(13,4.0);
    f->SetParameter(14,0.0018315);
    f->SetParameter(15,0.00790974);
    f->SetParameter(16,1.56139);
    f->SetParameter(17,1.39);
    f->SetParameter(18,0.428328);
    f->SetParameter(19,0.00115601);
    f->SetParameter(20,-0.0849149);
    f->SetParameter(21,2.58396);
    f->SetParameter(22,0.419646);
    f->SetParameter(23,10.0);
    f->SetParameter(24,0.0124098);
    f->SetParameter(25,0.0780432);
    f->SetParameter(26,0.182612);
    f->SetParameter(27,0.144939);
    f->SetParameter(28,-0.435278);
    f->SetParameter(29,0.0642399);
    f->SetParLimits(0, -10, 10);
    f->SetParLimits(1, -10, 10);
    f->SetParLimits(2,-10,10);
    f->SetParLimits(3,-10,10);
    f->SetParLimits(4,-10,10);
    f->SetParLimits(5, -10, 10);
    f->SetParLimits(6, -10, 10);
    f->SetParLimits(7,-10,10);
    f->SetParLimits(8,-10,10);
    f->SetParLimits(9,-10,10);
    f->SetParLimits(10, -10, 10);
    f->SetParLimits(11, -10, 10);
    f->SetParLimits(12,-10,10);
    f->SetParLimits(13,-10,10);
    f->SetParLimits(14,-10,10);
    f->SetParLimits(15, -10, 10);
    f->SetParLimits(16, -10, 10);
    f->SetParLimits(17,-10,10);
    f->SetParLimits(18,-10,10);
    f->SetParLimits(19,-10,10);
    f->SetParLimits(20, -10, 10);
    f->SetParLimits(21, -10, 10);
    f->SetParLimits(22,-10,10);
    f->SetParLimits(23,-10,10);
    f->SetParLimits(24,-10,10);
    f->SetParLimits(25, -10, 10);
    f->SetParLimits(26, -10, 10);
    f->SetParLimits(27,-10,10);
    f->SetParLimits(28,-10,10);
    f->SetParLimits(29,-10,10);*/
    
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
    for (int i=0; i < step_distance; i++){
        for (int j=0; j < step_angle; j++){
            if (count[i][j]<=0||isnan(Energy[i][j]/count[i][j])||abs(Energy[i][j])>100) {
                h2D->Fill(dis[i],ang[j],rr);
                x << dis[i] << endl;
                y << ang[j] << endl;
                z << rr << endl;
                xyz << dis[i] << " " << ang[j] << " " << rr << endl;
                continue;
            }
                h2D->Fill(dis[i],ang[j],Energy[i][j]/count[i][j]);
                x << dis[i] << endl;
                y << ang[j] << endl;
                z << Energy[i][j]/count[i][j] << endl;
                xyz << dis[i] << " " << ang[j] << " " << Energy[i][j]/count[i][j] << endl;
                //xyz << "(" << dis[i] << "," << ang[j] << "," << Energy[i][j]/count[i][j] << "),";
                rr = Energy[i][j]/count[i][j];
        }
    }
    
    cout << "Passed:" << cunt << " / " << cunt0 << endl;
    c1->cd();
    h2D->Draw("CONT4Z");
    //h2D->Draw("surf3");
    //h2D->Fit(f,"R");
    //f->Draw("SAME");
    file.close();
    x.close();
    y.close();
    z.close();
    xyz.close();
    return 0;
}
