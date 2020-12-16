Double_t par[8] = {1.18772, 2.58453, 0.138447, 1.01062, 0.0261287, 0.845129, 0.389123, 0.688391};

Double_t FFF(const Double_t *x){
    Double_t r = sqrt(x[0]*x[0]+x[1]*x[1]+0.00001);
    return (exp(-par[1]*r)+par[2]*exp(-par[3]*r)+par[4]*exp(-par[5]*r))/pow(r,par[7]);
    //return exp(-1/cos(x[0]))+0*x[1];
}
Double_t MMM(Double_t distance, Double_t angle){
    //time_t begin,end;
    //begin = clock();
    if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
    else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    angle *= TMath::DegToRad();
    Double_t L0 = par[0];
    Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    double a[2] = {xm,ym};
    double b[2] = {xp,yp};
    //const double ERRORLIMIT = 1E+2;
    ROOT::Math::Functor wf(&FFF,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE,0.0001,0.0001,1000);
    ig.SetFunction(wf);
    Double_t value = ig.Integral(a,b);
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return par[6]*value;
}

struct CZ{
    std::vector<pair<Double_t, Double_t>> pars;
    const Double_t distance_step = 0.1;
    CZ(){
        std::string in_name1 = "doc/p1.txt";
        std::ifstream p1;
        p1.open(in_name1,std::ios::in);
        std::string in_name2 = "doc/p2.txt";
        std::ifstream p2;
        p2.open(in_name2,std::ios::in);
        std::string str1,str2;
        while (std::getline(p1, str1)&&std::getline(p2, str2)) {
            std::stringstream strStream1(str1);
            std::stringstream strStream2(str2);
            float par1, par2, nl;
            strStream1 >> nl >> par1;
            strStream2 >> nl >> par2;
            //cout << par1 << ", " << par2 << endl;
            pars.push_back(make_pair(par1,par2));
        }
        p1.close();
        p2.close();
    }
    Double_t newfunc(Double_t distance, Double_t angle){
        Double_t value1, value2;
        if (distance < 0.5*distance_step){
            Double_t a = (0.5*distance_step - distance)/distance_step;
            value1 = pars[0].first;//-a*(pars[1].first-pars[0].first);
            value2 = pars[0].second;//-a*(pars[1].second-pars[0].second);
        }else if (distance >= 0.5*distance_step && distance < 3.5){
            int index = (int)((distance-0.5*distance_step)/distance_step);
            Double_t a = ((distance-0.5*distance_step) - index*distance_step)/distance_step;
            value1 = pars[index].first+a*(pars[index+1].first-pars[index].first);
            value2 = pars[index].second+a*(pars[index+1].second-pars[index].second);
        }else {
            return exp(-1.25*distance);
        }
        return value2*angle+value1;
    }
}method;

int test2(){
    int bin1(100),bin2(100),bin3(100);
    float tx(800),ty(600);
    //double xmin(0),xmax(3.5),ymin(-0.6),ymax(0.6),zmin(0),zmax(0.1);
    double xmin(0),xmax(3.5),ymin(0),ymax(1),zmin(0),zmax(90);
    
    TCanvas* c1=new TCanvas("PANDA1","test1",tx,ty);
    TCanvas* c2=new TCanvas("PANDA2","test2",tx,ty);
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
    
    TH1D* h1D = new TH1D("Hist1","h1",30,0,300);
    h1D->SetLineColor(kBlack);
    h1D->SetLineWidth(2);
    h1D->GetXaxis()->SetTitle("Time");
    h1D->GetYaxis()->SetTitle("Enties");
    //h1D->GetZaxis()->SetTitle("E");
    h1D->GetXaxis()->CenterTitle();
    h1D->GetYaxis()->CenterTitle();
    
    TH3D* h3D = new TH3D("Hist2","h2",bin1,xmin,xmax, bin3,zmin,zmax, bin2,ymin,ymax);
    h3D->SetMarkerStyle(7);
    h3D->SetMarkerColorAlpha(kAzure+3, 0.5);
    //h2D->GetYaxis()->SetTitle("E_{ci}");h2D->GetXaxis()->SetTitle("E_{truth}");
    h3D->GetZaxis()->SetTitle("E_{digi}");h3D->GetXaxis()->SetTitle("distance");h3D->GetYaxis()->SetTitle("angle");
    //h1D->GetZaxis()->SetTitle("E");
    h3D->GetXaxis()->CenterTitle();
    h3D->GetYaxis()->CenterTitle();
    h3D->GetZaxis()->CenterTitle();
    
    TGraph2D *g = new TGraph2D();
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetZaxis()->SetTitle("E_{digi}");
    g->GetXaxis()->SetTitle("distance");
    g->GetYaxis()->SetTitle("angle");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->GetZaxis()->CenterTitle();
    
    time_t begin,end,t(0.0);
    int N(0);
    Double_t x(0), y(0);
    for (int i = 0; i < 100; i++){
        cout << i << "%" << endl;
        y = 0;
        for (int j = 0; j < 100; j++){
            begin = clock();
            //if (Value(x,y)>1) cout << "XXX" << endl;
            //Double_t value = Value(x,y);
            //Double_t value = VALUE(x,y);
            //Double_t value = MYVALUE(x,y);
            //Double_t value = result(x,y);
            //Double_t value = SHOWER_DIGI(x,y);
            //Double_t value = Shower_Function.shower_Digi(x,y);
            //Double_t value = MMM(x,y);
            Double_t value = method.newfunc(x,y);
            //Double_t value = result_test(x);
            end = clock();
            h1D->Fill(end-begin);
            //if ( value>1 || value<0 ) cout << "ERROR!!" << endl;
            t += (end - begin);
            g->SetPoint(N,x,y,value);
            y += 0.45;
            N++;
        }
        x += 0.035;
    }
    cout << "TIME:" << t/CLK_TCK << endl;
    
    /*TGraph *g = new TGraph();
    int N(0);
    for (double p = 0.3; p < 3; p+=0.01){
        for(int i=0;i<1600;i++){
            g->SetPoint(N,i*0.001, mf(i*0.001,p));
            //g->SetPoint(N,i*0.001, ff(i*0.001,p));
            N++;
        }
    }*/
    
    
    TF1 *f=new TF1("f","[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[2]*x*x",0.001,0.75);
    c2->cd();
    g->Draw("p.");
    c1->cd();
    h1D->Draw();
    //g->Fit(f,"R");
    //h3D->Draw("SCAT");
    
    return 0;
}
