#ifndef KNN1_POLYNOMREGRESSION_H
#define KNN1_POLYNOMREGRESSION_H

#include <iostream>
#include "matrix.h"

using namespace std;

class PolynomRegression {
public:

    PolynomRegression(unsigned int mA);
    void setXandT(std::vector<double>xA,std::vector<double>tA);
    void computeAandBandW(void);
    inline double y(double xA);
    double error(void);


private:

    std::vector<double> xInE;
    std::vector<double> tInE;
    knn::matrix AE;
    knn::matrix bE;
    knn::matrix wE;
    unsigned int mE;
};

PolynomRegression::PolynomRegression(unsigned int mA)
{
  mE = mA;
}


void PolynomRegression::setXandT(std::vector<double> xA, std::vector<double> tA)
{
    xInE=xA;
    tInE=tA;
}

double PolynomRegression::y(double xA)
{
    double y = 0;
//aequivalent zu dem Summenzeichen die for-Schleife
    for (unsigned int m = 0; m <= mE; m++)
       y += wE(m+1,1)*pow(xA,m); // Matrix soll mit 1,1 beginnen und wE hat nur die erste Spalte
    return y;
}

double PolynomRegression::error(void)
{
    double errorsum = 0;

//Iteration über alle p-Werte
    for(int index = 0; index != tInE.size(); index++)
    errorsum += pow(y(xInE[index])- tInE[index],2); 
// xInE[index] steht für xp

    return errorsum;
}

void PolynomRegression::computeAandBandW(void)
{  //A
    AE = knn::matrix(mE+1,mE+1,0);

         vector<double>::iterator index;

   for (unsigned int i = 1; i <= mE+1; i++)
   {
       for (unsigned int j = 1; j <= mE+1; j++)
       {
            for(index = xInE.begin(); index != xInE.end(); index++)
               AE(i,j) += pow(*index,(j+i-2));
       }
   }
  //b
  bE =  knn::matrix(mE+1, 1, 0);
  for (unsigned int i = 1; i <= mE+1; i++)
   {

            for(unsigned int n = 0; n != xInE.size() ; n++)
               bE(i,1) += tInE[n] * pow(xInE[n],i-1);
       }
//Die zuerst invertierte Matrix AE mit bE multiplizieren um w-Werte zu erhalten
    AE.invert();
    wE = AE*bE; // muss wE vorher initialisiert werden ?, eigentlich ist AE*bE bereits eindeutig definiert... 




}

#endif KNN1_POLYNOMREGRESSION_H
