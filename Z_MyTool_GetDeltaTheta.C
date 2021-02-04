string Z_MyTool_GetDeltaTheta(double Theta_cent = 181, double Phi_Range = 9){
    double Phi_cent = 0.5;
    Phi_Range = fabs(Phi_Range);
    if (Theta_cent < 0 || fabs(Theta_cent)>180.001){
        cout << "Input theta center:(DEG)" << endl;
        cin >> Theta_cent;
        cout << "Input Phi range:(DEG)   default: 9 DEG" << endl;
        cin >> Phi_Range;
        if (Theta_cent <= 0) {cout << "-E  ERROR!!" << endl;return "XXXXXXXXX";}
        if (Phi_Range <= 0) Phi_Range = 9;
    }
    cout << "-INFO  Theta Center:  " << Theta_cent << " DEG" << endl;
    cout << "-INFO  Phi Center:  " << Phi_cent << " DEG" << endl;
    cout << "-INFO  Phi Range:  " << Phi_Range << " DEG" << endl;
    double Theta_cent_save = Theta_cent;
    double Phi_Range_save = Phi_Range;
    if (Theta_cent == 90) Theta_cent = 89.99999;
    Theta_cent *= TMath::DegToRad();
    Phi_Range *= TMath::DegToRad();
    TF1 *f = new TF1("func","tan([0]+x)-tan([0]-x)",0,Phi_Range);
    f->SetParameter(0, fabs(90*TMath::DegToRad() - Theta_cent));
    double DeltaTheta = f->GetX(Phi_Range);
    DeltaTheta *= TMath::RadToDeg();
    cout << "-INFO  Result => Theta Range: " << Theta_cent_save << " +- " << DeltaTheta << "  (" << 2*DeltaTheta << ")" << endl;
    cout << "******************" << endl << endl;
    cout << "tht(" << Theta_cent_save-DeltaTheta << ", " << Theta_cent_save+DeltaTheta << ")" << ":"
         << "phi(" << Phi_cent-Phi_Range_save/2 << ", " << Phi_cent+Phi_Range_save/2 << ")" << endl << endl;
    cout << "******************" << endl;
    double theta_m = Theta_cent_save-DeltaTheta, theta_p = Theta_cent_save+DeltaTheta;
    double phi_m = Phi_cent-Phi_Range_save/2, phi_p = Phi_cent+Phi_Range_save/2;
    ostringstream out1,out2,out3,out4;
    out1<<theta_m;
    out2<<theta_p;
    out3<<phi_m;
    out4<<phi_p;
    string str_theta_m = out1.str(), str_theta_p = out2.str(), str_phi_m = out3.str(), str_phi_p = out4.str();
    return "tht(" + str_theta_m + ", " + str_theta_p + ")" + ":"
    + "phi(" + str_phi_m + ", " + str_phi_p + ")" ;
}
