int Read_v1(){
    int bin1(100),bin2(200),bin3(50);
    float tx(800),ty(600);
    double xmin(0),xmax(3),ymin(0),ymax(45),zmin(0),zmax(1.01);
    const int step_distance(100),step_angle(100);
    
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
    
    TH3D* h2D = new TH3D("hvx0vy0","vx vs vy",bin1,xmin,xmax,bin2,ymin,2*ymax,bin3,zmin,zmax);
    //TH2D* h2D = new TH2D("hvx0vy0","vx vs vy",200,0,5,200,0,1);
    h2D->GetYaxis()->SetTitle("angle");
    h2D->GetXaxis()->SetTitle("d (cm)");
    h2D->GetZaxis()->SetTitle("E");
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    h2D->GetZaxis()->CenterTitle();
    h2D->SetMarkerStyle(7);
    h2D->SetMarkerColorAlpha(kAzure+3, 1);
    
    TF2 *f=new TF2("f","[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*y*y*y+[6]*y*y*y*y+[7]*y*y*y*y*y",0,3);
    //f->SetLineWidth(2);
    //f->SetLineColor(kRed);
    f->SetParameters(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
    //f->SetParLimits(0, -100, 100);
    //f->SetParLimits(1, -100, 100);
    //->SetParameters(2,-10,10);
    //f->SetParameters(3,-10,10);
    //f->SetParameters(4,-10,10);
    
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
    
    for (int i=0; i < step_distance; i++){
        for (int j=0; j < step_angle; j++){
            if (count[i][j]<=0||isnan(Energy[i][j]/count[i][j])||abs(Energy[i][j])>100) continue;
            if (j%2 == 0){
                h2D->Fill(dis[i],ang[j],Energy[i][j]/count[i][j]);
                x << dis[i] << endl;
                y << ang[j] << endl;
                z << Energy[i][j]/count[i][j] << endl;
                xyz << dis[i] << " " << ang[j] << " " << Energy[i][j]/count[i][j] << endl;
                //xyz << "(" << dis[i] << "," << ang[j] << "," << Energy[i][j]/count[i][j] << "),";
            }
            else{
                h2D->Fill(dis[i],90-ang[j],Energy[i][j]/count[i][j]);
                x << dis[i] << endl;
                y << 90-ang[j] << endl;
                z << Energy[i][j]/count[i][j] << endl;
                xyz << dis[i] << " " << 90-ang[j] << " " << Energy[i][j]/count[i][j] << endl;
                //xyz << "(" << dis[i] << "," << 90-ang[j] << "," << Energy[i][j]/count[i][j] << "),";
            }
        }
    }
    
    cout << "Passed:" << cunt << " / " << cunt0 << endl;
    c1->cd();
    h2D->Draw("SCAT");
    //h2D->Fit(f,"R");
    //f->Draw("SAME");
    file.close();
    x.close();
    y.close();
    z.close();
    xyz.close();
    return 0;
}
