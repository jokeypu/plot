int test1(){
TVector3 distance,pos(0,1,0),mom(0,0,1),vz(0,0,1);
                distance.SetPtThetaPhi(pos.Mag()*cos(mom.Angle(pos)), mom.Theta(), mom.Phi());
                distance = pos - distance;
		Double_t angle = distance.Angle(mom);
                //distance.SetPtThetaPhi(1,1.57,1.57);
		cout <<"Pt:"<< distance.Pt() << "Theta:"<<distance.Theta()<<"Phi:"<<distance.Phi()<< endl;
		cout <<"X:"<< distance.X() << "Y:"<<distance.Y()<<"Z:"<<distance.Z()<< endl;
		cout <<57.3* angle << endl;
return 0;
}
