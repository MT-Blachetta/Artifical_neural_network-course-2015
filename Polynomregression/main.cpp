#include "matrix.h"
#include "PolynomRegression.h"
#include "math.h"
#define _USE_MATH_DEFINES


// von Michael Blachetta, Martikelnummer: 108011206680


/* Antwort zu 1.1c 
1. Von Modellkomplexität 0 bis 10 sinkt der Restfehler stetig wo er bei M = 10 sein Minimum erreicht => offenbar Idealer M-Wert für die Funktion
2.Bei den Trainingsbeispielen sinkt der Restfehler viel schneller bis dahin als auf der Testmenge was jedoch zu erwarten ist da das System nicht darauf trainiert wurde
3. Für M > 10 scheint der Wert unerklärlicherweise direkt, massiv und überproportional anzusteigen. Auf der Testmenge noch viel stärker. Allerdings habe ich bei einigen Kommilitonen beobachtet das sie etwa ab M=9 unterschiedliche Werte haben. Vielleicht hat eine Variable einen Überlauf oder es wurde anders gerundet.
4. Etwa ab M = 13 scheint der Restfehler der Trainingsmenge zum selben Niveau zu konvergieren wie M <= 2, also wie zu Beginn */

/* Antwort zu 1.1e
1. Für die Trainingsmenge ist ein abweichendes jedoch ähnliches Verhalten zu beobachten wie bei der Sinusfunktion,
jedoch sind die Abfälle und Steigungen nicht so kontinuierlich und ruckartiger, das Minimum ist ebenfalls bei M = 10
2. Die Werte für die Testbeispiele scheinen bis M = 7 zu passen, dann aber steigen sie unmittelbar massiv an was 
unrealistisch ist, trotzdem habe ich mir viel Mühe gegeben es richtig zu programmieren und weiss nicht wie es dazu kommt */  
 

//zu Aufgabe 1.1 c

//Berechnet die Werte xp der Trainingsdaten
std::vector<std::vector<double> > computeTrainingSin(unsigned int pA)
{
    vector<vector<double> > trainingdata;

    double xp = 0;
    double tp = 0;

    vector<double> x;
    vector<double> t;

    for(unsigned int p = 1 ; p <= pA ; p++)
    {
        xp = ((2*p-1)*M_PI)/pA - M_PI + 0.01;
        x.push_back(xp);
        tp = sin(xp/2);
        t.push_back(tp);

    }

    trainingdata.push_back(x);
    trainingdata.push_back(t);

    return trainingdata;

}

//Berechnet die Werte tp der Testdaten
std::vector<std::vector<double> > computeTestSin(unsigned int pA)
{
    vector<vector<double> > testdata;

    double xp = 0;
    double tp = 0;

    vector<double> x;
    vector<double> t;

//beginnt mit 1 wie in Aufgabe
    for(unsigned int p = 1 ; p <= pA ; p++)
    {
        xp = (2*p*M_PI)/pA - M_PI + 0.01;
        x.push_back(xp);
        tp = sin(xp/2);
        t.push_back(tp);
    }

    testdata.push_back(x);
    testdata.push_back(t);

    return testdata;

}


//zu Aufgabe 1.1 e
std::vector<std::vector<double> > computeTrainingSinc(unsigned int pA)
{
    vector<vector<double> > trainingdata;

    double xp = 0;
    double tp = 0;

    vector<double> x;
    vector<double> t;

    for(unsigned int p = 1 ; p <= pA ; p++)
    {
        xp = ((2*p-1)*M_PI)/pA - M_PI + 0.01;
        x.push_back(xp);
        tp = sin(4*xp)/xp;
        t.push_back(tp);

    }

    trainingdata.push_back(x);
    trainingdata.push_back(t);

    return trainingdata;
}

std::vector<std::vector<double> > computeTestSinc(unsigned int pA)
{
    vector<vector<double> > testdata;

    double xp = 0;
    double tp = 0;

    vector<double> x;
    vector<double> t;

    for(unsigned int p = 1 ; p <= pA ; p++)
    {
        xp = (2*p*M_PI)/pA - M_PI + 0.01;
        x.push_back(xp);
        tp = sin(xp*4)/xp;
        t.push_back(tp);
    }

    testdata.push_back(x);
    testdata.push_back(t);

    return testdata;
}


int main(int argc, char** argv) {
    knn::init();

    //number of training examples: parameter 1, if given, else 11
    unsigned int P = argc > 1 ? (unsigned)atoi(argv[1]) : 11;
    //maximum polynomal degree: parameter 2, if given, else 19
    unsigned int MMax = argc > 2 ? (unsigned)atoi(argv[2]) : 19;

    std::ofstream sinOutL;
    sinOutL.open("SinError.txt",std::ios::out);
    std::ofstream sincOutL;
    sincOutL.open("SincError.txt",std::ios::out);


    for (unsigned int mL=0;mL<=MMax;++mL)
    {
        PolynomRegression regressorL(mL);
        std::vector<std::vector<double> > trainingDataL=computeTrainingSin(P);
        regressorL.setXandT(trainingDataL[0],trainingDataL[1]);
        regressorL.computeAandBandW();
        sinOutL<<mL<<"	"<<regressorL.error()<<"	";
        std::vector<std::vector<double> > testDataL=computeTestSin(P);
        regressorL.setXandT(testDataL[0],testDataL[1]);
        sinOutL<<regressorL.error()<<std::endl;
    }

    for (unsigned int mL=0;mL<=MMax;++mL)
    {
        PolynomRegression regressorL(mL);
        std::vector<std::vector<double> > trainingDataL=computeTrainingSinc(P);
        regressorL.setXandT(trainingDataL[0],trainingDataL[1]);
        regressorL.computeAandBandW();
        sincOutL<<mL<<"	"<<regressorL.error()<<"	";
        std::vector<std::vector<double> > testDataL=computeTestSinc(P);
        regressorL.setXandT(testDataL[0],testDataL[1]);
        sincOutL<<regressorL.error()<<std::endl;
    }
    sinOutL.close();
    sincOutL.close();

    return 0;
}
