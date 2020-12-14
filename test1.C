Double_t par(0);

void SetPar(Double_t pp){
    par = pp;
}
Double_t f(const Double_t *x){
    //return exp(-1*c[0]*d/c[3])+c[1]*exp(-1*c[2]*d/c[3]);
    return exp(-1*par/sin(x[0]))+0*x[1];
    //return exp(-1/cos(x[0]))+0*x[1];
}
Double_t ff(const Double_t x, const Double_t p){
    //return exp(-1*c[0]*d/c[3])+c[1]*exp(-1*c[2]*d/c[3]);
    return exp(-1*p/sin(x));
    //return exp(-1/cos(x[0]))+0*x[1];
}
Double_t mf(Double_t distance, Double_t p){
    time_t begin,end;
    begin = clock();
    SetPar(p);
    double a[2] = {0,0};
    double b[2] = {distance,1};
    const double ERRORLIMIT = 1E-0;
    ROOT::Math::Functor wf(&f,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    Double_t value = ig.Integral(a,b);
    end = clock();
    cout << "TIME:" << end - begin << endl;
    return value;
}

Double_t re_t(int n,double a,double b){
    double a_n[2] = {a,0};
    double b_n[2] = {b,1};
    const double ERRORLIMIT = 1E-3;
    ROOT::Math::Functor wf(&f,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    return ig.Integral(a_n,b_n);
}

Double_t re(int n,double a,double b){
    double sum1(0),sum2(0),x(a);
    double step = (b-a)/(2*n);
    for (int i = 1; i < 2*n; i+=2) sum1 += exp(-1*par/sin(x+i*step));
    for (int i = 2; i < 2*n-1; i+=2) sum2 += exp(-1*par/sin(x+i*step));
    return ((b-a)/(6*n))*(exp(-1*par/sin(a))+exp(-1*par/sin(b))+4*sum1+2*sum2);
}

/*Double_t myfunc_o(Double_t d,Double_t a,Double_t fPar,Double_t L0){
    time_t begin,end;
    begin = clock();
    if ( a > 90 && a <= 180 ) a = 180 - a;
    if ( a > 45 && a <= 90 ) a = 90 - a;
    a *= TMath::DegToRad();
    //double Lx(L0),Ly(L0);
    Double_t xp1,xm1,xp2,xm2;
    Double_t x0 = fabs(d*cos(a)), y0 = fabs(d*sin(a));
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    TVector2 point0(xp,ym), point1(xm,ym),  point2(xp,yp), point3(xm,yp);
    Double_t a0 = point0.Phi(), a1 = point1.Phi(), a2 = point2.Phi(), a3 = point3.Phi();
    if ((xp*xm<0) && (yp*ym<0)) {
        //ym = fabs(ym);
        //xm = fabs(xm);
        //cout << xp << ", " << xm << ", " << yp << ", " << ym << endl;
        //cout << a0 << ", " << a1 << ", " << a2 << ", " << a3 << endl;
        cout << Int(TMath::TwoPi()+TMath::PiOver2()-a0,TMath::PiOver2()+a2,fPar*xp) << endl;
        Double_t value = Int(TMath::PiOver2()+a1,TMath::PiOver2()+a3,fPar*xm)-Int(TMath::TwoPi()+TMath::PiOver2()-a0,TMath::PiOver2()+a2,fPar*xp)-Int(a2,a3,fPar*yp)+Int(a0,a1,fPar*ym);
        end = clock();
        cout << "TIME:" << end - begin << endl;
        return (TMath::TwoPi()+value)/fPar;
    }else{
        a0 = TMath::ATan(ym/xp);a1 = TMath::ATan(ym/xm);a2 = TMath::ATan(yp/xp);a3 = TMath::ATan(yp/xm);
        Double_t value = Int(TMath::PiOver2()-a0,TMath::PiOver2()-a2,fPar*xp)-Int(TMath::PiOver2()-a1,TMath::PiOver2()-a3,fPar*xm)-Int(a2,a3,fPar*yp)+Int(a0,a1,fPar*ym);
        end = clock();
        cout << "TIME:" << end - begin << endl;
        return value/fPar;
    }
}*/
/*Double_t deltaphi = fabs(point1.DeltaPhi(point2));
Double_t alpha_min = a1<=a2 ? a1 : a2;
Double_t result1 = deltaphi - Int(alpha_min+TMath::PiOver2(), alpha_min+deltaphi+TMath::PiOver2(), fPar*xp);
deltaphi = fabs(point2.DeltaPhi(point3));
alpha_min = a2<=a3 ? a2 : a3;
Double_t result2 = deltaphi - Int(alpha_min, alpha_min+deltaphi, fPar*yp);
deltaphi = fabs(point3.DeltaPhi(point4));
alpha_min = a3<=a4 ? a3 : a4;
cout << "a3:" << a3 << ", a4:" << a4 << endl;
cout << alpha_min+TMath::PiOver2() << "::" << alpha_min+deltaphi+TMath::PiOver2() << endl;
if (alpha_min+deltaphi+TMath::PiOver2() > TMath::TwoPi())
Double_t result3 = deltaphi - Int(alpha_min+TMath::PiOver2(), alpha_min+deltaphi+TMath::PiOver2(), fPar*xm);
cout << "re:" << result3 << endl;
deltaphi = fabs(point4.DeltaPhi(point1));
alpha_min = a4<=a1 ? a4 : a1;
Double_t result4 = deltaphi - Int(alpha_min, alpha_min+deltaphi, fPar*ym);*/

/*Double_t Int(Double_t a,Double_t b,Double_t p){
     time_t begin,end;
     begin = clock();
     //int n = 3;
     //Double_t ka = (a/TMath::PiOver2())-round(a/TMath::PiOver2());
     //Double_t kb = (b/TMath::PiOver2())-round(b/TMath::PiOver2());
     //Double_t sum1(0),sum2(0),x(a);
     //Double_t step = (b-a)/(2*n);
     //for (int i = 1; i < 2*n; i+=2) sum1 += exp(-1*p/sin(x+i*step));
     //for (int i = 2; i < 2*n-1; i+=2) sum2 += exp(-1*p/sin(x+i*step));
     //return ((b-a)/(6*n))*(exp(-1*p/sin(a))+exp(-1*p/sin(b))+4*sum1+2*sum2);
    if (sin((a+b)/2) < 0 ) {
        Double_t c = fabs(a-TMath::TwoPi());
        a = fabs(b-TMath::TwoPi());
        b = c;
        p = fabs(p);
    }
    Int_t n = 3;
    Double_t x;
    Double_t value1(0), value2(0);
    //0~0.72*sqrt(xm) n=60;
    Double_t step;
    //if ( p< -1) cout << "p<-1 !!!" << endl;
    if ( fabs(p) < 0.03 ){
        Int_t N = 20;
        Double_t a_cut = 0.3367*pow(p,0.4);
        x = a;
        step = (a_cut-a)/N;
        for (int i = 1; i < N; i++) value1 += exp(-p/sin(x+i*step));
        value1 = step*(exp(-p/sin(a))+exp(-p/sin(a_cut))/2+value1);
        a = a_cut;
        a_cut = 2.59111*pow(p,0.46);
        x=a_cut;
        step = (b-a_cut)/N;
        for (int i = 1; i < N; i++) value2 += exp(-p/sin(x+i*step));
        value2 = step*(exp(-p/sin(a_cut))+exp(-p/sin(b))/2+value2);
        b = a_cut;
    }
    x = a;
    Double_t sum1(0),sum2(0);
    step = (b-a)/(2*n);
    for (int i = 1; i < 2*n; i+=2) sum1 += exp(-p/sin(x+i*step));
    for (int i = 2; i < 2*n-1; i+=2) sum2 += exp(-p/sin(x+i*step));
    //if (fabs(a-xm)<0.01) return ((b-a)/6/n)*(TMath::ACos(xm/b)*exp(-ci*b)+4*sum1+2*sum2);else
    double vv = value1+value2+((b-a)/(6*n))*(exp(-p/sin(a))+exp(-p/sin(b))+4*sum1+2*sum2);
    end = clock();
    cout << "TIME:" << end - begin << endl;
    return vv;
}*/

Double_t ppp(0);
void Setppp(Double_t pp){
    ppp = pp;
}
Double_t fff(const Double_t *x){
    return exp(-ppp/sin(x[0]))+0*x[1];
}
Double_t Int(double a,double b,double p){
    //time_t begin,end;
    //begin = clock();
    Setppp(p);
    double a_n[2] = {a,0};
    double b_n[2] = {b,1};
    const double ERRORLIMIT = 1E-3;
    ROOT::Math::Functor wf(&fff,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    double vv = ig.Integral(a_n,b_n);
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return vv;
}

Double_t myfunc(Double_t d,Double_t a,Double_t fPar,Double_t L0){
    if ( a > 90 && a <= 180 ) a = 180 - a;
    if ( a > 45 && a <= 90 ) a = 90 - a;
    a *= TMath::DegToRad();
    Double_t x0 = d*cos(a), y0 = d*sin(a);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    TVector2 point1(xp,ym), point2(xp,yp),  point3(xm,yp), point4(xm,ym);
    Double_t a1 = point1.Phi(), a2 = point2.Phi(), a3 = point3.Phi(), a4 = point4.Phi();
    int w3(-1),w4(-1);
    if (xm>=0 && ym<0) w4 = 1;
    else if (xm<0 && ym<=0) w3 = w4 = 1;
    Double_t result1, result2, result3, result4;
    if(ym<0) {
        w4 = 1;
        result1 = a2 - (a1-TMath::TwoPi()) - Int(a1-TMath::TwoPi()+TMath::PiOver2(), a2+TMath::PiOver2(), fPar*xp);
        if (xm>0) result3 = a3 - (a4-TMath::TwoPi()) - Int(a4-TMath::TwoPi()+TMath::PiOver2(), a3+TMath::PiOver2(), fPar*xm);
        else {
            w3 = 1;
            result3 = a4 - a3 - Int(a3+TMath::PiOver2(), a4+TMath::PiOver2(), fPar*xm);
        }
        result2 = a3 - a2 - Int(a2, a3, fPar*yp);
        result4 = a1 - a4 - Int(a4, a1, fPar*ym);
    }else{
        result1 = a2 - a1 - Int(a1+TMath::PiOver2(), a2+TMath::PiOver2(), fPar*xp);
        result2 = a3 - a2 - Int(a2, a3, fPar*yp);
        result3 = a3 - a4 - Int(a4+TMath::PiOver2(), a3+TMath::PiOver2(), fPar*xm);
        result4 = a4 - a1 - Int(a1, a4, fPar*ym);
    }
    
    
    
    //if (fabs(a3 - a2) > TMath::Pi()) a3<a2 ? a2 -= TMath::TwoPi() : a2 += TMath::TwoPi();
    //Double_t result2 = a3 - a2 - Int(a2, a3, fPar*yp);
    //if (fabs(a4 - a3) > TMath::Pi()) a4<a3 ? a3 -= TMath::TwoPi() : a3 += TMath::TwoPi();
    //Double_t result3 = a4 - a3 - Int(a3+TMath::PiOver2(), a4+TMath::PiOver2(), fPar*xm);
    //if (ym>0){
       // Double_t result1 = a2 - a1 - Int(a1+TMath::PiOver2(), a2+TMath::PiOver2(), fPar*xp);
       // Double_t result4 = a1 - a4 - Int(a4, a1, fPar*ym);
    //}
    //if (fabs(fabs(a4 - a3)-TMath::Pi()) < 0.1) result3=0.0;
    //else
    //result3 = fabs(xm)*(a3 - a4 - Int(-TMath::ATan(L0/(2*0.15))+TMath::PiOver2(), TMath::ATan(L0/(2*0.15))+TMath::PiOver2(), fPar*0.15))/0.15;
    //if (a4 < a3) result3 = a3 - a4 - Int(a4+TMath::PiOver2(), a3+TMath::PiOver2(), fPar*fabs(xm));
    //else result3 = a4 - a3 - Int(a3+TMath::PiOver2(), a4+TMath::PiOver2(), fPar*fabs(xm));
    //if (fabs(a1 - a4) > TMath::Pi()) a1<a4 ? a4 -= TMath::TwoPi() : a4 += TMath::TwoPi();
    //Double_t result4;
    //if (a1 < a4 ) result4 = a4 - a1 - Int(a1, a4, fPar*ym);
    //else result4 = a1 - a4 - Int(a4, a1, fPar*ym);
    return (result1+result2+w3*result3+w4*result4)/fPar;
}

Double_t result(Double_t distance,Double_t angle){
    //time_t begin,end;
    //begin = clock();
    Double_t p[7] = {1.22069, 37147,  168.255,  -15302.4, 44.5161, 1489.1, 1.74012};
    Double_t value =(p[1]*myfunc(distance,angle,p[2],p[0])+p[3]*myfunc(distance,angle,p[4],p[0])+p[5]*myfunc(distance,angle,p[6],p[0]))/(3816*0.2*TMath::TwoPi());
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return value;
}

Double_t result_test(Double_t distance){
    //time_t begin,end;
    //begin = clock();
    Double_t value = exp(-1.25*distance);
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return value;
}
/*Double_t NFunc(Double_t d,Double_t a,Double_t fPar,Double_t L0){
    if ( a > 90 && a <= 180 ) a = 180 - a;
    if ( a > 45 && a <= 90 ) a = 90 - a;
    a *= TMath::DegToRad();
    Double_t x0 = d*cos(a), y0 = d*sin(a);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    Double_t min_distance = xm, max_distance = sqrt(xp*xp+yp*yp);
    if (xm>0 && ym>0) min_distance = sqrt(xm*xm+ym*ym);
    else if (xm<0 && ym<0) {
        return 0;
        //min_distance = sqrt(xm*xm+ym*ym);
        //return ((L0+xm)*L0/(TMath::Pi()*max_distance*max_distance))*(TMath::TwoPi()*(1-exp(-fPar*max_distance)))/fPar+(-xm*L0/(TMath::Pi()*min_distance*min_distance))*TMath::TwoPi()*(1-exp(-fPar*min_distance))/fPar;
        
    }
    Double_t sum =  TMath::TwoPi()*(exp(-fPar*min_distance)-exp(-fPar*max_distance))/fPar;
    return sum*L0*L0/(TMath::Pi()*(max_distance*max_distance-min_distance*min_distance));
}*/

Double_t NFunc(Double_t x,Double_t y,Double_t fPar,Double_t L0){
    Double_t x0, y0;
    if (fabs(x)<fabs(y)) {x0 = fabs(y);y0 = fabs(x);}
    else {x0 = fabs(x);y0 = fabs(y);}
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    Double_t min_distance = xm, max_distance = sqrt(xp*xp+yp*yp);
    if (xm>0 && ym>0) min_distance = sqrt(xm*xm+ym*ym);
    else if (xm<0) {
        return 2*L0*L0*(1-exp(-fPar*max_distance))/(max_distance*max_distance*fPar);
        //return TMath::TwoPi()*(1-exp(-fPar*max_distance))/fPar;
        //min_distance = fabs(xm);
        //return TMath::TwoPi()*(1-exp(-fPar*min_distance)+(exp(-fPar*min_distance)-exp(-fPar*max_distance)) * (L0*L0/TMath::Pi()-min_distance*min_distance)/(max_distance*max_distance-min_distance*min_distance))/fPar;
    }
    return 2*L0*L0*(exp(-fPar*min_distance)-exp(-fPar*max_distance))/((max_distance+min_distance)*(max_distance-min_distance)*fPar);
}

Double_t Exec(Double_t distance, Double_t angle, Double_t fPar, Double_t L0){
    Int_t N = 2;
    if (distance < 3*L0){
        if (distance < 1.42*L0) N = 18;
        else N = 10;
    }
    //if (distance <= L0*sqrt(2)) N = (int)25*exp(-0.3*fabs(distance));
    //if (distance <= L0/2) N = 25;
    Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    //N = (int)50*exp(-(1.2-angle/80)*fabs(x0));
    //else N = (int)30*exp(-(1.2-angle/80)*fabs(x0));
    //if (distance < 3*L0  &&  distance > L0*sqrt(2)) N = (int)8*exp(-0.8*distance);
    //if (distance <= L0*sqrt(2)) N = (int)50*exp(-0.8*fabs(xm));
    Double_t step = 2*L0/N, step_over2 = L0/N, sum = 0;
    Double_t xi = xm+step_over2;
    for (int i = 0; i < N; i++){
        Double_t yi = ym+step_over2;
        for (int j = 0; j < N; j++){
            sum += NFunc(xi,yi,fPar,step_over2);
            yi += step;
        }
        xi += step;
    }
    return sum;
}

Double_t Value(Double_t distance,Double_t angle){
    //time_t begin,end;
    //begin = clock();
    if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
    if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    angle *= TMath::DegToRad();
    Double_t p[7] = {1.22069, 1405,  3.15,  169, 0.887, 45, 0.354};
    Double_t value =(p[1]*Exec(distance,angle,p[2],p[0])+p[3]*Exec(distance,angle,p[4],p[0])+p[5]*Exec(distance,angle,p[6],p[0]))/(3816*0.2*TMath::TwoPi());
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return value;
}

Double_t NEWF(Double_t x, Double_t y, Double_t par){
    Double_t  r = sqrt(x*x+y*y)+0.1;
    return exp(-par*r)/r;
}
Double_t INT(Double_t x0, Double_t L0, Double_t y, Double_t par){
    //Double_t L0 = (b-a)/2;
    Double_t a = x0 - L0, b = x0 + L0;
    if (a*b<0){
        a *= -1;
        if (b < a) {
            Double_t ys = b;
            b = a;
            a = ys;
        }
        Double_t s = b-a;
        Double_t result1 = a/18*(NEWF(0,y,par)+NEWF(a,y,par)+4*(NEWF(a/6,y,par)+NEWF(a/2,y,par)+NEWF(5*a/6,y,par))+2*(NEWF(a/3,y,par)+NEWF(2*a/3,y,par)));
        Double_t result2 = s/18*(NEWF(a,y,par)+NEWF(b,y,par)+4*(NEWF(a+s/6,y,par)+NEWF(a+s/2,y,par)+NEWF(a+5*s/6,y,par))+2*(NEWF(a+s/3,y,par)+NEWF(a+2*s/3,y,par)));
        return 2*result1+result2;
    }
    return L0/9*(NEWF(a,y,par)+NEWF(b,y,par)+4*(NEWF(a+L0/3,y,par)+NEWF(a+L0,y,par)+NEWF(a+5*L0/3,y,par))+2*(NEWF(a+2*L0/3,y,par)+NEWF(a+4*L0/3,y,par)));
}
Double_t RES(Double_t distance, Double_t angle, Double_t par, Double_t L0){
    Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    if (ym*yp<0){
            //2*(0~ym)+(ym~yp)
            ym *= -1;
            if (yp < ym) {
                Double_t ys = yp;
                yp = ym;
                ym = ys;
            }
            Double_t s = yp-ym;
            Double_t result1 = ym/18*(INT(x0,L0,0,par)+INT(x0,L0,ym,par)+4*(INT(x0,L0,ym/6,par)+INT(x0,L0,ym/2,par)+INT(x0,L0,5*ym/6,par))+2*(INT(x0,L0,ym/3,par)+INT(x0,L0,2*ym/3,par)));
            Double_t result2 = s/18*(INT(x0,L0,ym,par)+INT(x0,L0,yp,par)+4*(INT(x0,L0,ym+s/6,par)+INT(x0,L0,ym+s/2,par)+INT(x0,L0,ym+5*s/6,par))+2*(INT(x0,L0,ym+s/3,par)+INT(x0,L0,ym+2*s/3,par)));
            return 2*result1+result2;
        
    }
    return L0/9*(INT(x0,L0,y0-L0,par)+INT(x0,L0,y0+L0,par)+4*(INT(x0,L0,y0-2*L0/3,par)+INT(x0,L0,y0,par)+INT(x0,L0,y0+2*L0/3,par))+2*(INT(x0,L0,y0-L0/3,par)+INT(x0,L0,y0+L0/3,par)));
}
Double_t VALUE(Double_t distance,Double_t angle){
    //time_t begin,end;
    //begin = clock();
    if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
    if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    angle *= TMath::DegToRad();
    Double_t p[7] = {1.22069, 1405,  3.15,  169, 0.887, 45, 0.354};
    Double_t value =(p[1]*RES(distance,angle,p[2],p[0])+p[3]*RES(distance,angle,p[4],p[0])+p[5]*RES(distance,angle,p[6],p[0]))/(3816*0.2*TMath::TwoPi());
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return value;
}

Double_t DIGITAL_ARC(Double_t a, Double_t b, Double_t xm, Double_t ci){
    Int_t n = 4;
    Double_t x;
    Double_t value1(0),value2(0);
    //0~0.72*sqrt(xm) n=60;
    Double_t step;
    if ( xm < 0.1 && b < 2.5 ){
        Int_t N = 10;
        Double_t a_cut = 0.72*sqrt(xm);
        x = a;
        step = (a_cut-a)/N;
        for (int i = 1; i < N; i++) value1 += TMath::ACos(xm/(x+i*step))*exp(-ci*(x+i*step));
        value1 = step*((TMath::ACos(xm/a)*exp(-ci*a)+TMath::ACos(xm/a_cut)*exp(-ci*a_cut))/2+value1);
        a = a_cut;
        
        /*a_cut = 1.35*sqrt(xm);;
        x = a_cut;
        step = (b-a_cut)/N;
        for (int i = 1; i < N; i++) value2 += TMath::ACos(xm/(x+i*step))*exp(-ci*(x+i*step));
        value2 = step*((TMath::ACos(xm/a_cut)*exp(-ci*a_cut)+TMath::ACos(xm/b)*exp(-ci*b))/2+value2);
        b = a_cut;*/
        
        /*Int_t N = 20;
        Double_t a_cut = 0.3367*pow(p,0.4);
        x = a;
        step = (a_cut-a)/N;
        for (int i = 1; i < N; i++) value1 += exp(-p/sin(x+i*step));
        value1 = step*(exp(-p/sin(a))+exp(-p/sin(a_cut))/2+value1);
        a = a_cut;
        a_cut = 2.59111*pow(p,0.46);
        x=a_cut;
        step = (b-a_cut)/N;
        for (int i = 1; i < N; i++) value2 += exp(-p/sin(x+i*step));
        value2 = step*(exp(-p/sin(a_cut))+exp(-p/sin(b))/2+value2);
        b = a_cut;*/
    }
    
    x = a;
    Double_t sum1(0),sum2(0);
    step = (b-a)/(2*n);
    for (int i = 1; i < 2*n; i+=2) sum1 += TMath::ACos(xm/(x+i*step))*exp(-ci*(x+i*step));
    for (int i = 2; i < 2*n-1; i+=2) sum2 += TMath::ACos(xm/(x+i*step))*exp(-ci*(x+i*step));
    //if (fabs(a-xm)<0.01) return ((b-a)/6/n)*(TMath::ACos(xm/b)*exp(-ci*b)+4*sum1+2*sum2);else
            return value1+value2+((b-a)/6/n)*(TMath::ACos(xm/a)*exp(-ci*a)+TMath::ACos(xm/b)*exp(-ci*b)+4*sum1+2*sum2);
}
Double_t SRT(Double_t distance, Double_t angle, Double_t ci, Double_t L0){
    int N(0);
    Double_t S[4];
    int w0(-1),w2(-1);
    Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
    for (int j = 0; j < 2; j++){
        Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
        if (j==0) {
            if (xm>=0 && ym<0) w2 = 1;
            else if (xm<0 && ym<=0) w0 = w2 = 1;
        }
        Double_t x_zd = fabs(xm);
        for (int i=0; i < 2; i++){
            Double_t r_min = sqrt(x_zd*x_zd+ym*ym), r_max = sqrt(x_zd*x_zd+yp*yp);
            Double_t theta_min= TMath::ACos(x_zd/r_min)  ,theta_max= TMath::ACos(x_zd/r_max);
            //if (r_min < 0.01) theta_min = 0;
            if (ym*yp<0)
                S[N] = theta_min/ci*(1-exp(-ci*r_min))+theta_max/ci*(1-exp(-ci*r_max))-2*DIGITAL_ARC(x_zd,r_min,x_zd,ci)-DIGITAL_ARC(r_min,r_max,x_zd,ci);
            else S[N] = (theta_max-theta_min+theta_min*exp(-ci*r_min)-theta_max*exp(-ci*r_max))/ci-DIGITAL_ARC(r_min,r_max,x_zd,ci);
            x_zd = xp;
            N++;
        }
        x0 = distance*sin(angle), y0 = distance*cos(angle);
    }
    return w0*S[0]+S[1]+w2*S[2]+S[3];
}

Double_t MYVALUE(Double_t distance,Double_t angle){
    //time_t begin,end;
    //begin = clock();
    if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
    if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    angle *= TMath::DegToRad();
    Double_t p[7] = {1.22069, 1405,  3.15,  169, 0.887, 45, 0.354};
    Double_t value =(p[1]*SRT(distance,angle,p[2],p[0])+p[3]*SRT(distance,angle,p[4],p[0])+p[5]*SRT(distance,angle,p[6],p[0]))/(3816*0.2*TMath::TwoPi());
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return value;
}

Double_t PPS[2];
void SetPPS(Double_t p0, Double_t p1){
    PPS[0] = p0;
    PPS[1] = p1;
}
Double_t FFUNC(const Double_t *x){
    return PPS[0]*exp(-PPS[1]*sqrt(x[0]*x[0]+PPS[0]*PPS[0]))/(PPS[0]*PPS[0]+x[0]*x[0]);
}
/*Double_t FUNC_INT(Double_t a, Double_t b, Double_t y, Double_t ci){
    SetPPS(y,ci);
    double a_n[2] = {a,0};
    double b_n[2] = {b,1};
    const double ERRORLIMIT = 1E+3;
    ROOT::Math::Functor wf(&FFUNC,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    return ig.Integral(a_n,b_n);
}*/
Double_t hs(Double_t x, Double_t y, Double_t ci){
    Double_t xi = x*x+y*y;
    return y*(1-exp(-ci*sqrt(xi)))/xi;
}
Double_t FUNC_INT(Double_t a, Double_t b, Double_t y, Double_t ci){
    Int_t N, n = 3;
    float value(0.0);
    //if (y<0.001) y = 0.001;
    if (y < 0.1){
        //cout << b-a << endl;
        N = 3+100*(b-a);
        float step = (b-a)/N;
        float k = a+step, m = b-step/2;
        //for (int i = 1; i < N; i++) value += hs(a+i*step,y,ci);
        for (float i = k; i < m; i+=step) value += hs(i,y,ci);
        value = step*((hs(a,y,ci)+hs(b,y,ci))/2+value);
    }else{
        float sum1(0),sum2(0);
        float step = (b-a)/(2*n);
        float Twostep = 2*step, StepOver2 = step/2;
        float k1 = a + step, k2 = k1 + step, m1 = b-StepOver2, m2 = m1 - StepOver2;
        for (float i = k1; i < m1; i+=Twostep) sum1 += hs(i,y,ci);
        for (float i = k2; i < m2; i+=Twostep) sum2 += hs(i,y,ci);
        //for (int i = 1; i < 2*n; i+=2) sum1 += hs(a+i*step,y,ci);
        //for (int i = 2; i < 2*n-1; i+=2) sum2 += hs(a+i*step,y,ci);
        value = step/3*(hs(a,y,ci)+hs(b,y,ci)+4*sum1+2*sum2);
    }
    //return (exp(-ci*a)-exp(-ci*b))/(ci*y)-value;
    return TMath::ATan(b/y)-TMath::ATan(a/y)-value;
}
/*Double_t FUNC_INT(Double_t a, Double_t b, Double_t y, Double_t ci){
    Int_t n = 10;
    Double_t x;
    Double_t value1(0),value2(0);
    //if (a==0) a+=0.00001;
    //0~0.72*sqrt(xm) n=60;
    Double_t step;
    if ( y < 0.1) n = 10000;
    if ( y < 0.1 && ci<-3000){
        Int_t N = 1000;
        Double_t a_cut = 1.6*y;
        x = a;
        step = (a_cut-a)/N;
        for (int i = 1; i < N; i++) value1 += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
        value1 = step*((y*exp(-ci*sqrt(a*a+y*y))/(a*a+y*y)+y*exp(-ci*sqrt(a_cut*a_cut+y*y))/(a_cut*a_cut+y*y))/2+value1);
        a = a_cut;
        
        a_cut = 6*y;
        x = a_cut;
        step = (b-a_cut)/N;
        for (int i = 1; i < N; i++) value2 += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
        value2 = step*((y*exp(-ci*sqrt(b*b+y*y))/(a*a+y*y)+y*exp(-ci*sqrt(a_cut*a_cut+y*y))/(a_cut*a_cut+y*y))/2+value2);
        b = a_cut;
    }
    x = a;
    Double_t sum1(0),sum2(0);
    step = (b-a)/(2*n);
    for (int i = 1; i < 2*n; i+=2) sum1 += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
    for (int i = 2; i < 2*n-1; i+=2) sum2 += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
    Double_t value = value1+value2+((b-a)/6/n)*(y*exp(-ci*sqrt(a*a+y*y))/(a*a+y*y)+y*exp(-ci*sqrt(b*b+y*y))/(b*b+y*y)+4*sum1+2*sum2);
    if (fabs(value) > 4)
    cout << value << "," << a << ", " << sum1 << ", " << sum2 << ", " << value1 << ", " << y << endl;
    //if (fabs(a-xm)<0.01) return ((b-a)/6/n)*(TMath::ACos(xm/b)*exp(-ci*b)+4*sum1+2*sum2);else
    return value;
    //int n = 3000;
    //double sum1(0),sum2(0),x(a);
    //double step = (b-a)/(2*n);
    //for (int i = 1; i < 2*n; i+=2) sum1 += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
    //for (int i = 2; i < 2*n-1; i+=2) sum2 += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
    //return ((b-a)/(6*n))*(y*exp(-ci*sqrt(a*a+y*y))/(a*a+y*y)+y*exp(-ci*sqrt(b*b+y*y))/(b*b+y*y)+4*sum1+2*sum2);
}*/
/*Double_t FUNC_INT(Double_t a, Double_t b, Double_t y, Double_t ci){
    Int_t N = 4000;
    Double_t x(a), value(0);
    Double_t step = (b-a)/N;
    for (int i = 1; i < N; i++) value += y*exp(-ci*sqrt((x+i*step)*(x+i*step)+y*y))/((x+i*step)*(x+i*step)+y*y);
    return step*((y*exp(-ci*sqrt(a*a+y*y))/(a*a+y*y)+y*exp(-ci*sqrt(b*b+y*y))/(b*b+y*y))/2+value);
}*/

/*Double_t FUNC_INT(Double_t a, Double_t b, Double_t y, Double_t ci){
    int n = 50;
    Double_t I(0), n1, n2;
    Double_t f_a = y*exp(-ci*sqrt(a*a+y*y))/(a*a+y*y);
    Double_t f_b = y*exp(-ci*sqrt(b*b+y*y))/(b*b+y*y);
    for (int i = 0; i < n; i++) {
        n1 = 0.01*(2*i+rand()/(RAND_MAX+1.));
        n2 = 0.01*(2*i+1+rand()/(RAND_MAX+1.));
        n1 = (b-a)*n1+a;
        n2 = (b-a)*n2+a;
        //cout << n1 << ", " << n2 << "   " << "(" << a << ", " << b << ")"<< endl;
        //(((b-a)*n1+a)-f_b)/(f_a-f_b);
        //((y*exp(-ci*sqrt(n2*n2+y*y))/(n2*n2+y*y)-f_b)/(f_a-f_b);
        I += (y*exp(-ci*sqrt(n2*n2+y*y))/(n2*n2+y*y) + y*exp(-ci*sqrt(n1*n1+y*y))/(n1*n1+y*y))-2*f_b;
    }
    I = 0.01*I;
    //if (f_a > f_b)
    return (b-a)*I + (b-a)*f_b;
    //else return (b-a)*(I + f_a);
}*/
Double_t LINE_INT(Double_t x_n, Double_t y_n, Double_t ci, Double_t L0){
    TVector2 p1(x_n-L0,y_n),p2(x_n+L0,y_n);
    Double_t value(0);
    if ( x_n-L0 > 0) value = FUNC_INT(x_n-L0,x_n+L0,y_n,ci);
    else value = 2*FUNC_INT(0,L0-x_n,y_n,ci)+FUNC_INT(L0-x_n,x_n+L0,y_n,ci);
    return (p1.Phi()-p2.Phi()-value)/ci;
}
Double_t BOLOCK(Double_t distance, Double_t angle, Double_t ci, Double_t L0){
    Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    //Double_t line_1[] = {xm,y0}, line_2[] = {xp,y0}, line_3[] = {x0,ym}, line_4[] = {x0,yp};
    int w1(1),w3(1);
    if (xm>0) {
        w1 = -1;
        if (ym>0) w3 = -1;
    }
    return w1*LINE_INT(y0,fabs(xm),ci,L0)+LINE_INT(y0,xp,ci,L0)+w3*LINE_INT(x0,fabs(ym),ci,L0)+LINE_INT(x0,yp,ci,L0);
    
}
Double_t SHOWER_DIGI(Double_t distance,Double_t angle){
    //time_t begin,end;
    //begin = clock();
    if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
    if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    angle *= TMath::DegToRad();
    Double_t p[7] = {1.22069, 1405,  3.15,  169, 0.887, 45, 0.354};
    Double_t value =(p[1]*BOLOCK(distance,angle,p[2],p[0])+p[3]*BOLOCK(distance,angle,p[4],p[0])+p[5]*BOLOCK(distance,angle,p[6],p[0]))/(3816*0.2*TMath::TwoPi());
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return value;
}

struct INTEGRAL {
    Double_t y, ci, L0;
    
    void SetPar(Int_t index, Double_t par){
        if (index == 1) L0 = par;
        else if (index == 2) ci = par;
        else if (index == 3) y = par;
        else std::cout << "Parameter error!!" << std::endl;
    }
    
    /*Double_t func(Double_t x){
        Double_t xi = x*x+y*y;
        return y*(1-exp(-ci*sqrt(xi)))/xi;
    }*/
    
    Double_t func_Int(Double_t a, Double_t b){
        float value(0);
        if (y < 0.1){
            Int_t N = 3+200*(b-a);
            float step = (b-a)/N, k = a+step, m = b-step/2;
            for (float i = k; i < m; i+=step) value += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            return step*((y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value);
        }else{
            Int_t n = 4;
            float step = (b-a)/(2*n), Twostep = 2*step, StepOver2 = step/2;
            float sum1(0), sum2(0), k1 = a + step, k2 = k1 + step, m1 = b - StepOver2, m2 = m1 - step;
            for (float i = k1; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            for (float i = k2; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            return step/3*(y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
        }
        //return TMath::ATan(b/y)-TMath::ATan(a/y)-value;
    }
    
    Double_t line_Int(Double_t x_n, Double_t y_n){
        //TVector2 p1(x_n-L0,y_n),p2(x_n+L0,y_n);
        //Double_t value(0);
        SetPar(3,y_n);
        if ( x_n-L0 > 0) return func_Int(x_n-L0,x_n+L0);
        else return 2*func_Int(0,L0-x_n)+func_Int(L0-x_n,x_n+L0);
        //return p1.Phi()-p2.Phi()-value;
        //return value;
    }
    
    /*Double_t line_Int(Double_t x_n, Double_t y_n){
        SetPar(3,y_n);
        Double_t a = fabs(x_n-L0), b = x_n+L0;
        float value(0);
        if (y < 0.1){
            Int_t N = 3+100*(b-a);
            float step = (b-a)/N, k = a+step, m = b-step/2;
            for (float i = k; i < m; i+=step) value += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            value = step*((y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value);
        }else{
            Int_t n = 3;
            float step = (b-a)/(2*n), Twostep = 2*step, StepOver2 = step/2;
            float sum1(0), sum2(0), k1 = a + step, k2 = k1 + step, m1 = b - StepOver2, m2 = m1 - step;
            for (float i = k1; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            for (float i = k2; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            value = step/3*(y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
        }
        if ( x_n-L0 > 0) return value;
        else {
            float value_add(0);
            b = a;
            if (y < 0.1){
                Int_t N = 3+100*b;
                float step = b/N, m = b-step/2;
                for (float i = step; i < m; i+=step) value_add += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
                value_add = step*((y*(1-exp(-ci*sqrt(y*y)))/(y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value_add);
            }else{
                Int_t n = 3;
                float step = b/(2*n), Twostep = 2*step, StepOver2 = step/2;
                float sum1(0), sum2(0), m1 = b - StepOver2, m2 = m1 - step;
                for (float i = step; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
                for (float i = Twostep; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
                value_add = step/3*(y*(1-exp(-ci*sqrt(y*y)))/(y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
            }
            return value + 2*value_add;
        }
    }*/
    
      /*Double_t line_Int(Double_t x_n, Double_t y_n){
        float value(0), value_add(0);
        SetPar(3,y_n);
        Double_t a = fabs(x_n-L0), b = x_n+L0;
        if (y < 0.1){
            Int_t N = 3+100*(b-a);
            float step = (b-a)/N, k = a+step, m = b-step/2;
            for (float i = k; i < m; i+=step) value += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            value = step*((y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value);
            if ( x_n-L0 > 0) return value;
            else{
                b = a;
                N = 3+100*b;
                step = b/N; m = b-step/2;
                for (float i = step; i < m; i+=step) value_add += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
                value_add = step*((y*(1-exp(-ci*sqrt(y*y)))/(y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value_add);
                return value+2*value_add;
            }
        }else{
            Int_t n = 3;
            float step = (b-a)/(2*n), Twostep = 2*step, StepOver2 = step/2;
            float sum1(0), sum2(0), k1 = a + step, k2 = k1 + step, m1 = b - StepOver2, m2 = b - 3*StepOver2;
            for (float i = k1; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            for (float i = k2; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            value = step/3*(y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
            if ( x_n-L0 > 0) return value;
            else{
                b = a;
                step = b/(2*n); Twostep = 2*step; StepOver2 = step/2;
                sum1 = 0; sum2 = 0; m1 = b - StepOver2; m2 = b - 3*StepOver2;
                for (float i = step; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
                for (float i = Twostep; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
                value_add = step/3*(y*(1-exp(-ci*sqrt(y*y)))/(y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
                return value+2*value_add;
            }
        }
       }*/
    
    Double_t block_Int(Double_t distance, Double_t angle, Double_t ci){
        Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
        Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
        SetPar(2,ci);
        Int_t w1(1),w3(1);
        if (xm>0) {
            w1 = -1;
            if (ym>0) w3 = -1;
        }
        return w1*line_Int(y0,fabs(xm))+line_Int(y0,xp)+w3*line_Int(x0,fabs(ym))+line_Int(x0,yp);
    }
    
    /*Double_t block_Int(Double_t distance, Double_t angle, Double_t ci){
        Double_t value(0);
        Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
        Double_t xx[4] = {x0-L0, x0+L0 , y0-L0, y0+L0};
        int w[4] = {1, 1, 1, 1};
        SetPar(2,ci);
        if (xx[0]>0) { w[0] = -1; if (xx[2]>0) w[2] = -1; }
        int n(0);
        for (int i = 0; i < 2; i++){
            int index = (n+2)%4;
            Double_t xxm = xx[index], fxxm = -xxm;
            Double_t xxp = xx[index+1];
            for (int j = 0; j < 2; j++) {
                SetPar(3,fabs(xx[n]));
                if ( xxm > 0) value += w[n]*func_Int(xxm,xxp);
                else value += w[n]*(2*func_Int(0,fxxm)+func_Int(fxxm,xxp));
                n++;
            }
        }
        return value;
        
    }*/
    
    Double_t shower_Digi(Double_t distance,Double_t angle){
        //time_t begin,end;
        //begin = clock();
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
        angle *= TMath::DegToRad();
        //Double_t p[7] = {1.22069, 1405,  3.15,  169, 0.887, 45, 0.354};
        //Double_t p[7] = {1.22,  3.15,  169.0/1405.0, 0.887, 45.0/1405.0, 0.354 ,0.357472};
        Double_t p[7] = {0.22, 396.258, 3.35543, 47.0842, 0.703498, 0.112794, -0.326313};
        SetPar(1,p[0]);
        //Double_t value =p[6]*(block_Int(distance,angle,p[1])/p[1]+p[2]*block_Int(distance,angle,p[3])/p[3]+p[4]*block_Int(distance,angle,p[5])/p[5])/3;//(3816*0.2*TMath::TwoPi());
        Double_t value =(p[1]*block_Int(distance,angle,p[2])/p[2]+p[3]*block_Int(distance,angle,p[4])/p[4]+p[5]*block_Int(distance,angle,p[6])/p[6])/(3816*0.2*TMath::TwoPi());
        //end = clock();
        //cout << "TIME:" << end - begin << endl;
        return value;
    }
    
}Shower_Function;

int test1(){
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
    
    TH1D* h1D = new TH1D("Hist1","h1",30,0,30);
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
        y = 0;
        for (int j = 0; j < 1000; j++){
            begin = clock();
            //if (Value(x,y)>1) cout << "XXX" << endl;
            //Double_t value = Value(x,y);
            //Double_t value = VALUE(x,y);
            //Double_t value = MYVALUE(x,y);
            //Double_t value = result(x,y);
            //Double_t value = SHOWER_DIGI(x,y);
            Double_t value = Shower_Function.shower_Digi(x,y);
            //Double_t value = result_test(x);
            end = clock();
            h1D->Fill(end-begin);
            //if ( value>1 || value<0 ) cout << "ERROR!!" << endl;
            t += (end - begin);
            g->SetPoint(N,x,y,value);
            y += 0.09;
            N++;
        }
        x += 0.03;
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
