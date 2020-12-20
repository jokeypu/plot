#include <iostream>
#include <string>
using namespace std;
int File_Process(std::string dir_name){
    std::string in_name = "doc/"+ dir_name +".txt";
    std::string out_name = "doc/"+ dir_name +"_R.txt";
    
    std::ifstream File_in;
    File_in.open(in_name, std::ios::in);
    std::ofstream File_out;
    File_out.open(out_name,std::ios::out);
    if (!File_in.is_open()) {
        std::cout << "read node data failed" << std::endl; return 1;
    }
    if (!File_out.is_open()) {
        std::cout << "wirte node data failed" << std::endl; return 1;
    }
    
    std::string str;
    //std::getline(File_in, str);
    
    float distance_max(8), angle_max(45);
    const int bin_x(80), bin_y(90);
    int bin_table[bin_x][bin_y] = {0};
    float Energy_table[bin_x][bin_y] = {0};
    float step_x = distance_max/bin_x, step_y = angle_max/bin_y;
    while (std::getline(File_in, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        if ( distance >= distance_max ) continue;
        if ( (angle_max-45)<0.01 && angle > 45 && angle <= 90 ) angle = 90 - angle;
        int ibinx = (int) (distance/step_x);
        int ibiny = (int) (angle/step_y);
        Energy_table[ibinx][ibiny] += energy;
        bin_table[ibinx][ibiny]++;
        //cout << distance << ", " << angle << ", " << energy << endl;
    }
    File_in.close();
    
    for (int i = 0; i < bin_x; i++){
        for (int j = 0; j < bin_y; j++){
            if (bin_table[i][j] < 3) continue;
            float d = step_x/2 + i*step_x, a = step_y/2 + j*step_y, E = Energy_table[i][j]/bin_table[i][j];
            File_out << d << " " << a << " " << E << endl;
        }
    }
    File_out.close();
    
    //TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
    /*TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
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
    gStyle->SetTitleOffset(1.0,"xyz");*/
    
    /*TGraph2D *g1 = new TGraph2D("doc/DigiEnergy.txt","%lg %lg %lg");
    g1->SetMarkerStyle(7);
    g1->SetMarkerColorAlpha(kAzure+3, 0.5);
    g1->GetZaxis()->SetTitle("E_{digi}");
    g1->GetXaxis()->SetTitle("distance");
    g1->GetYaxis()->SetTitle("angle");
    g1->GetXaxis()->CenterTitle();
    g1->GetYaxis()->CenterTitle();
    g1->GetZaxis()->CenterTitle();*/
    
    /*TGraph2D *g2 = new TGraph2D(out_name,"%lg %lg %lg");
    g2->SetMarkerStyle(7);
    g2->SetMarkerColorAlpha(kAzure+3, 0.5);
    g2->GetZaxis()->SetTitle("E_{digi}");
    g2->GetXaxis()->SetTitle("distance");
    g2->GetYaxis()->SetTitle("angle");
    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();
    g2->GetZaxis()->CenterTitle();
    
    //c1->cd();
    //g1->Draw("p.");
    c2->cd();
    g2->Draw("tri2 p0");*/
    return 0;
}
