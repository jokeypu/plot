/*Double_t myfunc(Double_t x) {
 Double_t p1(-0.2262), p2(-0.1186), x0(4.399);
 Double_t p0 = p1*p1 -p2;
 Double_t value;
 if (x <= x0) value = p0/((x-p1) * (x-p1) - p2);
 else{
 Double_t fx0 = p0/((x0-p1) * (x0-p1) - p2);
 Double_t c = fx0 * (x0-p1) * (x0-x) / (p0*1.1512925);
 value = fx0 * pow(10, c);
 }
 return value;
 }*/
/*
 Double_t myfunc(Double_t x) {
 Double_t value = exp(-1.176 * x);
 return value;
 }
 */
/*
 Double_t myfunc(Double_t x) {
 if (x<1.4) return 0.99;
 else return exp(-2.89 * sqrt(x-1.4));
 }
 */

Double_t myfunc(Double_t x, Double_t p0, Double_t p1) {
    if (x<1.7) return 0.33;
    else return exp(-1 * p0 * sqrt(x-p1));
}

/*
 Double_t myfunc(Double_t x) {
 if (x<1.7) return 0.33;
 else return exp(-1.25*x);
 }
 */
