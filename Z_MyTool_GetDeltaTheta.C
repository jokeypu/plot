string Z_MyTool_GetDeltaTheta(double Theta_cent = 181, double Phi_Range = 9){
    double Phi_cent = 0;
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
    return "tht(" + Theta_cent_save-DeltaTheta + ", " + Theta_cent_save+DeltaTheta + ")" + ":"
    + "phi(" + Phi_cent-Phi_Range_save/2 + ", " + Phi_cent+Phi_Range_save/2 + ")" ;
}
