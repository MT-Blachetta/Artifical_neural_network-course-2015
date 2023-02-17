#include<iostream>
#include<string>
#include<sstream>
#include<vector>
#include<math.h>
#include<fstream>
#include<map>
#include"matrix.h"
//#include"gnuplot-iostream.h"


using namespace std;

namespace bl {

//Aktivierungsfunktionen


struct Layer{

unsigned int neuronNr;
bool IN;

knn::matrix weights;
knn::matrix changes;

knn::matrix neuron;

};

struct DataSet{

vector< vector<double> > X;
vector< vector<double> > Y;

};

double SignumFunktion(double argument) {

return 1/(1 + exp(-argument));

}

double AbleitungSignumFunktion(double argument) {

double n = 1/(1 + exp(-argument));
return n - (n*n);
}

double SinusFunktion(double arg, unsigned int k, unsigned int fnr){

	if(fnr%2)
	return cos(k*arg);
	else
	return sin(k*arg);
	}


double SinusAbleitung(double arg, unsigned int k, unsigned int fnr){

	if(fnr%2)
	return -k*sin(k*arg);
    else
	return k*cos(k*arg);

}

double Hauptfunktion(vector<double> inputvektor){

double insum = 0;

for(int d = 0; d < inputvektor.size(); d++)
insum += inputvektor[d];

cout << "insum = \n" << insum << endl;

double result = (double)(((long int)insum % 2) * 2) - 1.0;
cout << "result = \n" << result << endl;
return result;

}



vector<double> BinaryfromDecimal(unsigned int p, unsigned int in){

//dezimale Zahl p soll in einen Binärvektor umgewandelt werden
//z.B. p = 3 in 1100...(binärzahl)
assert(p>=0);

cout << "BinaryfromDecimal Übergabewert " << p << endl;
vector<double> ret;
double h = 0;
p = floor(p);
cout << "p = " << p << endl;

cout << "in = " << in << endl;

for(int d = 0; d < in; d++){
cout << "h(alt) = " << h << endl;
h = (int)floor(p/pow(2,d)) % 2 ;
cout << "h(neu) = " << h << endl;
//cout << "ret(alt) = " << ret[d] << endl;
ret.push_back(h);
cout << "ret(neu) = " << ret[d] << endl;
}

cout << endl;

return ret;

}

class net                                                    {

public:

bl::Layer* network;
double maximum;


//gnuplot pterminal;

DataSet Testdata;
DataSet Trainingdata;

int superior;

DataSet RandomFunction;



unsigned int layerSize;

/* net(vector<unsigned int> NxM);
net(vector<unsigned int> NxM, double tmax);

void MainMethod(double low, double top, double lernrate, const char* ResultPath);

void calcActivity(vector<double> input);
void backpropagate(vector<double> eingabe, vector<double> soll, double nu);
void learnFunction(vector <vector<double> > arg, vector <vector<double> > goal, double learnRate, double epoch);
double meanError(vector <vector<double> > ein, vector <vector<double> > aus);
bool isReady(double ac, int lnr);
void printOutput(void);
void NetzAusgabe(vector<double> x);
void NetzInfo(void);
void fuctionToFile(int range, const char* title);
bool layerFileInit(const char* fname, unsigned int l);
//void layerFileClose(unsigned int l);
void layerpwl(unsigned int l, unsigned int t);
bool init_plotfile(const char* filename);
void println_plotfile(double index, double value);
void learnFunction(vector <vector<double> > arg, vector <vector<double> > goal, double learnRate, double epoch, const char* outfile);
void randInitWeights(double down, double up);
void createTrainingdata(unsigned int nr);
void createTestdata(double down, double up, unsigned int dataSize);
void randomTestdata(double down, double up); */

//net(){}


net(vector<unsigned int> NxM){




//plotfile.open("main.txt");
//cout << "NETWORK CONSTRUCTER: " << endl;
layerSize = NxM.size();
int insize = NxM[0];
maximum = pow(2,insize) - 1;
//cout << "layerSize = " << layerSize << endl;
//cout << endl;
network = new Layer[layerSize];

//cout << "layer: 0 " << endl;
network[0].IN = true;
//cout << "is(inputlayer) = " << network[0].IN << endl;
network[0].neuronNr = NxM[0];
//cout << "INPUT neuronNr = " << network[0].neuronNr << endl;
network[0].neuron = knn::matrix(2, NxM[0], 1);
//cout << "neuron.colSize = " << network[0].neuron.colSize() << endl;
//cout << "neuron.rowSize = " << network[0].neuron.rowSize() << endl;
//cout << endl;
for(int l = 1; l < layerSize; l++){
//cout << "layer: " << l << endl;
network[l].IN = false;
network[l].neuronNr = NxM[l];
network[l].neuron = knn::matrix(2, NxM[l], 1);
network[l].weights = knn::matrix(NxM[l], NxM[l-1]+1, 1);
network[l].changes = knn::matrix(NxM[l], NxM[l-1]+1, 0);
//cout << "NxM["  << l  << "] = " << NxM[l] << endl;
//cout << "NxM["  << l-1  << "] = " << NxM[l-1] << endl;
//cout << "is(inputlayer) = " << network[l].IN << endl;
//cout << "#neurons  = " << network[l].neuron.colSize() << endl;
//cout << "weights.colSize = " << network[l].weights.colSize() << endl;
//cout << "weights.rowSize = " << network[l].weights.rowSize() << endl;
//cout << endl;
}

}

net(vector<unsigned int> NxM, double tmax, string suffix){

cout << "open plotfile...\n";


string title = "ConstructorINFO " + suffix;
ofstream plotfile (title);

if(!plotfile.is_open())
    cout << "OPENING ERROR \n";

plotfile << "1. OK \n";

maximum = tmax-1;
int maxIN = (int)ceil(log2(maximum));
plotfile << "maxIN 1 = " << maxIN << endl;
vector<double> temp = BinaryfromDecimal(tmax, maxIN);
maxIN = temp.size();
plotfile << "maxIN 2 = " << maxIN << endl;

for(int a = 0; a < maxIN; a++)
plotfile << "temp Vector " << a << " , VALUE = " << temp[a] << endl;

if(maxIN < NxM[0]){
plotfile << "NETWORK CONSTRUCTER: " << endl;
layerSize = NxM.size();
maximum = pow(2,NxM[0]) - 1;
plotfile << "4. OK \n";
plotfile << "maximum = " << maximum << endl;
plotfile << endl;
network = new Layer[layerSize];

plotfile << "layer: 0 " << endl;
network[0].IN = true;
plotfile << "is(inputlayer) = " << network[0].IN << endl;
network[0].neuronNr = NxM[0];
plotfile << "INPUT neuronNr = " << network[0].neuronNr << endl;
network[0].neuron = knn::matrix(2, NxM[0], 1);
plotfile << "neuron.colSize = " << network[0].neuron.colSize() << endl;
plotfile << "neuron.rowSize = " << network[0].neuron.rowSize() << endl;
plotfile << endl;
} else {
  plotfile << "NETWORK CONSTRUCTER: " << endl;
  plotfile << "layer: 0 " << endl;
    layerSize = NxM.size();
    network = new Layer[layerSize];
    network[0].IN = true;
plotfile << "is(inputlayer) = " << network[0].IN << endl;
network[0].neuronNr = maxIN;
NxM[0] = maxIN;
plotfile << "maxIN (changed) = " << maxIN << endl;
plotfile << "INPUT neuronNr = " << network[0].neuronNr << endl;
network[0].neuron = knn::matrix(2, maxIN, 1);

}

for(int l = 1; l < layerSize; l++){
plotfile << "layer: " << l << endl;
network[l].IN = false;
network[l].neuronNr = NxM[l];
network[l].neuron = knn::matrix(2, NxM[l], 1);
network[l].weights = knn::matrix(NxM[l], NxM[l-1]+1, 1);
network[l].changes = knn::matrix(NxM[l], NxM[l-1]+1, 0.11);
plotfile << "NxM["  << l  << "] = " << NxM[l] << endl;
plotfile << "NxM["  << l-1  << "] = " << NxM[l-1] << endl;
plotfile << "is(inputlayer) = " << network[l].IN << endl;
plotfile << "#neurons  = " << network[l].neuron.colSize() << endl;
plotfile << "weights.colSize = " << network[l].weights.colSize() << endl;
plotfile << "weights.rowSize = " << network[l].weights.rowSize() << endl;
plotfile << endl;
}
cout << "vor Testdaten erzeugung \n";
createTestdata(0,maximum,tmax);
cout << "Testdaten erzeugt \n";
createTrainingdata((int)maximum/2);
cout << "ERFOLGREICH ! \n\n\n";

plotfile.close();

//cout << endl;

                             }

~net(){

}

void MainMethod(double low, double top, double lernrate, double number, string ResultPath){

randInitWeights(low,top);

string errName;
string functionFile;
//funktion unklar

errName = "meanErrors_" + ResultPath;
functionFile = "START_NetworkFunction_" + ResultPath;

FunctionToFile((int)maximum, functionFile);
//plotFile(functionFile,"TEST","input","output",1,2,(int)maximum,2);

double vor;
vor = meanError(Testdata.X,Testdata.Y);

NetzInfo();
SimplelearnFunction(Trainingdata.X, Trainingdata.Y, lernrate, number, errName);

functionFile = "RESULT_NetworkFunction_" + ResultPath;
FunctionToFile((int)maximum, functionFile);

NetzInfo();

double nach;
nach = meanError(Testdata.X,Testdata.Y);

cout << "Fehler  VORHER: " << vor << endl;
cout << "Fehler NACHHER: " << nach << endl;

return;

}

void createTrainingdata(unsigned int nr){

unsigned int inputSize = network[0].neuronNr;
vector<double> IN;
knn::initRNG(0);
double generator = 0;


for(int p = 0; p < nr; p++){
generator = knn::randomFromInterval(0,maximum);
IN = BinaryfromDecimal(generator,inputSize);
Trainingdata.X.push_back(IN);
vector<double> tdv;
tdv.push_back(Hauptfunktion(IN));
Trainingdata.Y.push_back(tdv);
}



}

void createTestdata(double down, double up, unsigned int dataSize){

cout << "Aufruf OK \n";
cout << endl;
unsigned int inputSize = network[0].neuronNr;
cout << " inputsize " << inputSize << endl;
vector<double> tempIN;
//cout << "tempIN Erstellt: " << tempIN[0] << endl;

knn::initRNG(0);
double rand = 0;

cout << "Random Value Init OK \n";

 ////////////////////////////////////////////////////////////////////////////////////////////////

//cout << "Testdata Y initialized " << endl;

for(int p = 0; p < dataSize; p++){

rand = knn::randomFromInterval(down,up);
cout << "random: " << rand << endl;
tempIN = BinaryfromDecimal(rand,inputSize);
cout << "BinaryfromDecimal SUCCESS\n";
cout << "tempIN[0] = " << tempIN[0] << endl;


Testdata.X.push_back(tempIN);
cout << "Testdata push back\n";
vector<double> tdv;
tdv.push_back(Hauptfunktion(tempIN));
Testdata.Y.push_back(tdv);
cout << "Hauptfunktion SUCCESS\n";
}

cout << "TESTDATA SUCCESS \n";

}


void calcActivity(vector<double> input) {



assert(input.size() == network[0].neuronNr);
double netin = 0;
int nnr;


for(unsigned int i = 1; i <= network[0].neuronNr; i++)

network[0].neuron(1, i) = input[i-1];

for(unsigned int l = 1; l < layerSize; l++){

nnr = network[l-1].neuronNr;

for(unsigned int ln = 1; ln <= network[l].neuronNr; ln++){

netin += network[l].weights(ln, 1); 

for(unsigned int n = 1; n <= nnr; n++){

netin += network[l-1].neuron(1, n) * network[l].weights(ln, n+1); 

}


if(l == layerSize-1)
network[l].neuron(1, ln) =  netin;
else{
network[l].neuron(1, ln) = SignumFunktion(netin);
network[l].neuron(2, ln) = AbleitungSignumFunktion(netin);
}

netin = 0;

                                                }
        

                                           }

                                       


}

//outputs

void FunctionToFile(int range, string title){

ofstream datei (title);

if(!datei.is_open()){
    cout << "FILE CAN NOT BE OPENED !" << endl;
    return;
}

vector<double> x;

if(range >= maximum)
range = maximum - 1;



for(int pos = 0; pos <= range; pos++){
x = BinaryfromDecimal(pos,network[0].neuronNr);
calcActivity(x);
datei << pos << "  " << network[layerSize-1].neuron(1, 1) << endl;
}

datei.close();

}

void printOutput(void){
cout << endl;
cout << "RESULT: " << endl;
cout << endl;
for(int y = 1; y <= network[layerSize-1].neuron.colSize(); y++)
cout << "y" << y << "  = " << network[layerSize-1].neuron(1, y) << endl;
cout << endl;
}

void NetzAusgabe(vector<double> x){
calcActivity(x);
printOutput();
}

void NetzInfo(void){

//In Datei statt "cout" Ausgeben

cout << endl;
cout << "NETWORK INFO: \n" << endl;
cout << "Network layerSize: " << layerSize << endl;
cout << endl;

cout << "current Activity: " << endl;
printOutput();
cout << endl;
cout << "#inputs: " << network[0].neuron.colSize() << endl;
cout << endl;

for(int L = 1; L < layerSize; L++){
cout << endl;
cout << endl;
cout << "LAYER " << L << endl;
cout << endl;


cout << "#Neurons: " << network[L].neuronNr << endl;
cout << endl;
cout << "Activity \n" << endl;

for(int m = 1; m <= network[L].neuron.colSize(); m++)
cout << "z(" << L << ")" << m << " = " << network[L].neuron(1, m) << endl;

cout << endl;

cout << "WEIGHTS layer " << L << endl;

for(int k = 1; k <= network[L].weights.rowSize(); k++){
cout << "INcoming weights of Neuron " << k << " :" << endl;
for(int w = 1; w <= network[L].weights.colSize(); w++)
cout << "w" << k << "," << w-1 << " = " << network[L].weights(k, w) << endl;
}

}

}



bool isReady(double ac, int lnr)                           {
// Prüft,ob nach einem Lernverfahren (Epoche) sich die Gewichte kaum verändert haben und der Algorithmus fertig ist
assert(lnr > 0);

ac = ac*ac;

double sqchange;

for(unsigned int z = 1; z <= network[lnr].changes.rowSize() ; z++){

for(unsigned int s = 1; s <= network[lnr].changes.colSize() ; s++){
sqchange = network[lnr].changes(z, s);
sqchange = sqchange*sqchange;

if(sqchange > ac)
return false;

//sqchange = 0;

                                                       }

                                                       }

return true;
                                                           }


void backpropagate(vector<double> eingabe, vector<double> soll, double nu){

 unsigned int lastindex = network[layerSize-1].neuronNr;

knn::matrix lastErrors = knn::matrix(1, lastindex, 0);


calcActivity(eingabe);

for(int o = 0; o < soll.size(); o++){

lastErrors(1, o+1) = network[layerSize-1].neuron(1, o+1) - soll[o];

}

double zm = 0;
double fsig = 0;
int nnr;

			             for(int l = layerSize-1; l > 0; l--) {

nnr = network[l-1].neuronNr;

knn::matrix newErrors = knn::matrix(1, nnr, 0);

for(int ln = 1; ln <= nnr; ln++){

if(network[l-1].IN)// nicht für Inputschicht
break;


for(int s = 1; s <= lastErrors.colSize(); s++){
fsig += lastErrors(1, s) * network[l].weights(s, ln+1); // vom Bias Neuron wird kein Fehlersignal berechnet

}

zm = network[l-1].neuron(1, ln);

newErrors(1, ln) = network[l-1].neuron(2, ln) * fsig;

fsig = 0;


                                }



for(int r = 1; r <= network[l].neuronNr; r++){

network[l].changes(r, 1) = -nu * lastErrors(1, r);

network[l].weights(r, 1) += network[l].changes(r, 1);


for(int c = 2; c <= network[l].changes.colSize(); c++){
network[l].changes(r, c) = -nu * network[l-1].neuron(1, c-1) * lastErrors(1, r);

network[l].weights(r, c) += network[l].changes(r, c);


               }

                                             }


if(network[l-1].IN)
break;



lastErrors = newErrors;



                                                                                            }
          
}


double meanError(vector <vector<double> > ein, vector <vector<double> > aus) {
assert(ein.size() == aus.size());

double E = 0;
double meanE = 0;
double sumError = 0;

for(int p = 0; p < ein.size(); p++){

vector<double> soll = aus[p];
vector<double> in = ein[p];
assert(soll.size() == network[layerSize-1].neuron.colSize());

calcActivity(in);


for(int oj = 0; oj < soll.size(); oj++){

E = network[layerSize-1].neuron(1, oj+1) - soll[oj];

sumError += E;

}

meanE += (sumError*sumError);
sumError = 0;

}

double P = (double)aus.size();

return meanE/P ; // ((double)aus.size()) ;


}



void SimplelearnFunction(vector <vector<double> > arg, vector <vector<double> > goal, double learnRate, double epoch, string outfile){

assert(arg.size() == goal.size());

//Datei öffnen/erstellen
ofstream file (outfile);
if(!file.is_open())
    return;

double ce;
bool complete;

int hx;

for(int e = 0; e < epoch; e++){


for(int v = 0; v < arg.size(); v++)
backpropagate(arg[v],goal[v],learnRate);


ce = meanError(Testdata.X,Testdata.Y);

file << e << "  " << ce << endl;

for(int netL = 1; netL < layerSize; netL++){
complete = isReady(0.01,netL);
if(!complete)
break;
}

if(complete){

file.close();
return;
}

}



file.close();

}


void randInitWeights(double down, double up){

knn::initRNG(0);
double rand = 0;

int rows;
int cols;

for(int h = 1; h < layerSize; h++){

rows = network[h].weights.rowSize();
cols = network[h].weights.colSize();

for(int r = 1; r <= rows; r++){
for(int c = 1; c <= cols; c++){
rand = knn::randomFromInterval(down,up);
network[h].weights(r, c) = rand;
}

}

}
                                        }


/*
void plotFile(string fname, string title, string xlabel, string ylabel, unsigned int x, unsigned int y, int max_x, int max_y){

Gnuplot pterminal;

int xscale = (int)ceil(max_x/10.0);
int yscale = (int)ceil(max_y/10.0);
string s = "plot";
string suffix = ".jpeg";

pterminal << "reset\n";

pterminal << "set term jpeg size 800 600", endl;
pterminal << "set xlabel " << xlabel << endl;
pterminal << "set ylabel " << ylabel << endl;
pterminal << "set key left top\n";
pterminal << "set yzeroaxis\n";
pterminal << "set yrange [-2:" << max_y << "]\n";
pterminal << "set xrange [0:" << max_x << "]\n";
pterminal << "set xtics " << xscale << endl;
pterminal << "set ytics " << yscale << endl;
pterminal << "set label " << fname << endl;
pterminal << "set sample 1000\n";
pterminal << "set output " << s << title << suffix << endl;
pterminal << "plot " << fname << " using " << x << ":" << y << " with linespoints lc 2 lw 1 title '" << title << "'\n";

pterminal << "reset\n";


} */



void randomTestdata(double down, double up){


//Testdata vorher initialisieren
// need X/Y testdata.size

unsigned int inputSize = network[0].neuronNr;

vector<double> tempIN;

knn::initRNG(0);
double rand = 0;

for(int p = 0; p < Testdata.X.size(); p++){
rand = knn::randomFromInterval(down,up);
tempIN = bl::BinaryfromDecimal(rand,inputSize);
Testdata.X[p] = tempIN;
Testdata.Y[p][0] = Hauptfunktion(tempIN);
}

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*

void RANDOM_DATASET_GENERATOR(string newname, int defArea, int sup){

//size 10000

ofstream f_out ("generatedFunction.txt");

assert(file.is_open());


knn::initRNG(0);

double z = 0;


for(int i = 0; i < defArea; i++){
vector<double> inp;
vector<double> out;
z = i/100;
inp.push_back(z);
f_out << z << "  ";
RandomFunction.X.push_back(inp);
z = knn::randomFromInterval(0,sup);
f_out << z << endl;
out.push_back(z);
RandomFunction.Y.push_back(z);
}

f_out.close();



}

double RANDOM_FUNCTION(vector<double> arg){

	double ValueIndex = arg[0];
	ValueIndex *= 100;
	ValueIndex = floor(ValueIndex);

	unsigned int index = (int)ValueIndex;

    return RandomFunction.Y.[index][0];




}


CREATE_TRAININGSDATA(double xdef, unsigned int dataSize){

	knn::initRNG(0);

	double z = 0;

	for(int d = 0; d < dataSize; d++){
		z = knn::randomFromInterval(0,xdef-1);
		z = floor(ValueIndex);
		vector<double> vv;
		vv.push_back(z/100);
		Trainingdata.X.push_back(vv);
		vv = RandomFunction.Y.[(int)z];
		Trainingdata.Y.push_back(vv);

	}



}

CREATE_TESTDATA(double xdef, unsigned int dataSize){

	knn::initRNG(0);

	double z = 0;

	for(int d = 0; d < dataSize; d++){
		z = knn::randomFromInterval(0,xdef-1);
		z = floor(ValueIndex);
		vector<double> vv;
		vv.push_back(z/100);
		Testdata.X.push_back(vv);
		vv = RandomFunction.Y.[(int)z];
		Testdata.Y.push_back(vv);

	}



} */



};

}
