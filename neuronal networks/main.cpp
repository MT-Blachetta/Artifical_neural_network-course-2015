#include"net.h"
#include<vector>
#include<iostream>
#include<string>
#include<sstream>

/* Michael Blachetta 108011206680
   Negah Moradi 108012203826

Alle Aufgaben von Übungsblatt 1-5 können im Prinzip mit diesem Programm gelöst werden.

5.2

a) Implementieren sie das Zweischichtennetzwerk ->[im User interface können sie beliebig viele Layer samt hidden neurons (Modellkomplexität) erstellen, siehe "net" Konstruktor]
inklusive Error Backpropagation-Algorithmus ->[siehe "backpropagate(vector<double>,vector<double>,double)" und "SimplelearnFunction"] mit der in Aufgabe 5.1 berechneten Aktivierungsregeln ->[siehe "SignumFunktion(double)" und "backpropagation"]
mit variablem D ->[es kann mit "createTestdata" und in der UI eine int zahl angegeben werden die eine betrimmte Testdaten Menge erzeugt, M kann über Konstruktor variiert werden] und M

b) Erzeugen Sie alle moglichen Eingabebeispiele.->("BinaryfromDecimal" und "Hauptfunktion") Verwenden Sie zur Erzeugung der einzelnen Bei- ¨
spiele die in Aufgabe 5.1 c) entwickelte Berechnungsvorschrift. Teilen Sie die Menge der Beispiele
sinnvoll in Trainings- und Testmenge auf

c) Initialisieren Sie die Gewichte des Netzwerkes mit auf dem Intervall [−2, 2] ->("randomWeightInit") gleichverteilten Zufallszahlen.
Trainieren Sie das Netzwerk mit Backpropagation ->("simplelearnFunction"), einer quadratischen Fehlerfunktion ->("meanError")
und einer Lernrate η = 0.01 ->(jede beliebige Lernrate ist mögkich) in einem Online-Verfahren fur ca. 10000 Epochen. ->(wird im UI abgefragt als "learning cycles") Berechnen sie am ¨
Ende jeder Epoche den Fehler EMSSE auf der Trainings- und der Testmenge. ->(wird automatisch in Datei ausgegeben) Geben Sie diese
Fehler in eine Datei aus, sodass Sie den Fehlerverlauf plotten konnen. ->(ACHTUNG: das ist der eigentliche Grund für die Fehler, den Frühstadium des Programms und späte Abgabe da ich versucht habe zusätzliche Funktionalität einzufügen die Daten automatisch plottet, ich habe es bis jetzt hin so weit geschafft das ich eine "boost" libary installieren muss damit mein Code das könnte, da ich aber c++ und compiler Interna kaum kenne erweist es sich als schwierig)

d) Studieren Sie den Fehlerverlauf in Trainings- und Testmenge fur¨ D = 2, 5, 8, 10 und jeweils
mehreren Werten von M (auch großer als ¨ D). Fur welche der angegebenen Werte von ¨ D kann sich
das Netzwerk gut an die Funktion anpassen und warum? Wie hangt die Antwort von der Anzahl ¨
der verdeckten Neuronen ab?
->(das Programm ist dafür gemacht automatisch für beliebige Eingaben Netzfunktionen und Ausgabedaten zu erstellen und auch zu plotten, ohne manuellen Aufwand)



*/

int main(void){

unsigned int inputN = 2;
unsigned int outputN = 1;

unsigned int hidLayer = 0;
unsigned short hidbool = 1;

unsigned int NeuronNr = 2;

string fileSuffix;

double low = -1.0;
double top = 1.0;

unsigned int dataSize = 10;

vector<unsigned int> in;

cout << "Welcome to Neural Network Simulator 'PROTO' (alpha version) " << endl;
cout << endl;
cout << endl;

cout << "type in file output name suffix (string): \n";
cin >> fileSuffix;
cout << endl;

if(fileSuffix.empty())
fileSuffix = "DEFAULT.txt";

cout << endl;
cout << endl;

cout << "Choose Network Properties \n" << endl;
cout << endl;
cout << "How many inputs ? type an integer: \n";
cin >> inputN;
cout << endl;
assert(inputN > 0);
in.push_back(inputN);

cout << "Output number is here const set to 1 \n" << endl;


 /*cout << "How many network outputs ? type an integer: \n";
cin >> outputN;
cout << endl;
assert(outputN > 0); */

cout << "Add hidden Layers: \n\n";
while(hidbool){
cout << "#hidden Layer: " << hidLayer << endl;
cout << "for additional Layer type '1' for no type '0' " << endl;
cin >> hidbool;
cout << endl;

if(hidbool){
cout << "how many Neurons ? (integer): \n";
cin >> NeuronNr;
cout << endl;
in.push_back(NeuronNr);
}
else
cout << "Number hidden Layers: " << hidLayer << endl;

hidLayer++;

};

in.push_back(outputN);
cout << " \n" << endl;



cout << "Network Code \n";

for(int l = 0; l < in.size(); l++)
cout << in[l] << "  ";

cout << " \n" << endl;



cout << "Select Test and Training-data size; minimum 3 , type an integer: \n";
cin >> dataSize;
cout << endl;

if(dataSize < 3)
dataSize = 3;



cout << endl;
cout << "creating network... \n";
cout << endl;
cout << endl;


bl::net netz = bl::net(in,dataSize,fileSuffix); //--------------------------------------------------------------------

cout << "\n" << endl;
cout << "Random Initialisation; type in lower bound as x.x : \n";
cin >> low;
cout << endl;
cout << "Random Initialisation; type in upper bound as x.x : \n";
cin >> top;
cout << endl;
if(top <= low){
cout << "ERROR ! upper bound is to low \n";
return 0;}

cout << endl;
cout << "initializing weights random from " << low << " to " << top << endl;
cout << endl;
cout << endl;
cout << endl;


double nu = 0.1;

cout << "\nStarting learning process, type an learning rate: \n";
cin >> nu;
cout << endl;

unsigned int epoc;

cout << "Starting learning process, type number of learning cycle (int): \n";
cin >> epoc;
cout << endl;

cout << endl;
cout << endl;

cout << "Backpropagation START... \n";
cout << endl;
cout << endl;

netz.MainMethod(low,top,nu,epoc,fileSuffix);

//netz.MainMethod(-4,4,-0.01,100,"FehlerVerlauf.txt");

return 0;
//cout << " 1. Ok " << endl;

/*

vector<vector<double> > f;
vector<double> a = {0.0, 0.0};
vector<double> b = {1.0, 0.0};
vector<double> c = {0.0, 1.0};
vector<double> d = {1.0, 1.0};
f.push_back(a);
f.push_back(b);
f.push_back(c);
f.push_back(d);

for(int a = 0; a < f.size(); a++){
cout << " f vector Wert  " << a << endl;
for(int b = 0; b < f[a].size(); b++)
cout << " f[" << a << "]" << "[" << b << "] = " << f[a][b] << endl;
}

cout << endl;

vector< vector<double> > soll = {{-1.0},{1.0},{1.0},{-1.0}};


for(int a = 0; a < soll.size(); a++){
cout << " soll vector Wert  " << a << endl;
for(int b = 0; b < soll[a].size(); b++)
cout << " soll[" << a << "]" << "[" << b << "] = " << soll[a][b] << endl;
}

for(int k = 0; k < 4; k++)
netz.NetzAusgabe(f[k]);

//cout << " 3. Ok " << endl;
cout << endl;

netz.NetzInfo();

cout << endl;
cout << endl;

double t = 0;
t = netz.meanError(f,soll);

//cout << "Fehler VORHER: " << t << endl;



netz.learnFunction(f,soll,-0.01,100);

netz.NetzInfo();

double r = 0;
r = netz.meanError(f,soll);

cout << endl;

cout << endl;
cout << endl;

cout << "Fehler  VORHER: " << t << endl;
cout << endl;
cout << "Fehler NACHHER: " << r << endl;
//netz.test(1);

*/

}

/* das müsste eigentlich beliebig vieles von selbst plotten

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


}

*/
