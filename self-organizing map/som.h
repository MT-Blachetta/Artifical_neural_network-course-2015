#include<iostream>
#include<fstream>
#include<math.h>
#include"matrix.h"

using namespace std;

namespace som {

/* Eine 3-dimensionale Oberfläche soll mittels einer topologieerhaltenden 2-dimensionalen Kohonen-Karte
mit R x S Neuronen z(r,s) optimal abgedeckt werden. Jedes Neuron z(r,s) besitzt einen 3-dimensionalen
Zentrumsvektor µ(r,s)
, der die Position des Neurons im Eingaberaum beschreibt.


*/

struct matrixcluster{

unsigned int R; //es ist hier auch int möglich aber um einheitliche Typen zu verwenden ist es unsigned
unsigned int S;

/*
"soll mittels einer topologieerhaltenden 2-dimensionalen Kohonen-Karte
mit R x S Neuronen z(r,s) optimal abgedeckt werden" ->

jede Matrix hat einen Zeilenindex r = rows
und einen eindeutigen Spalteninde s = cols
R = matrix.getRowSize();
S = matrix.getColSize();

*/
knn::matrix first;
knn::matrix second;
knn::matrix third;

/*
"Jedes Neuron z(r,s) besitzt einen 3-dimensionalen
Zentrumsvektor µ(r,s), der die Position des Neurons im Eingaberaum beschreibt" ->

first = erste Komponente des Vektors µ_x1
second = zweite Komponente des Vektors µ_x2
third = dritte Komponente des Vektors µ_x3

*/


};


class neuroGrid {


public:

matrixcluster map;

double fnPara;

neuroGrid(double fnp){

fnPara = fnp;
/*
"Schreiben Sie eine Funktion, die eine Kohonen-Karte mit 1600 Neuronen zrs mit r; s = 0; ... ; 39
erzeugt und die Zentren µ(r,s)
r, s als zweidimensionales Gitter initialisiert" 
->
 */

double counter1 = 0;
double counter2 = 0;
map.R = 40;
map.S = 40;
map.first = knn::matrix(40, 40, 0);
map.second = knn::matrix(40, 40, 0);
map.third = knn::matrix(40, 40, 0); // <-  "während µ3 = 0 einen festen Wert bekommt"

/*

"wobei die Komponenten
µ1 und µ2 jeweils gleichmäßig im Bereich [-2; 2] liegen"

->

*/

for(unsigned int r = 1; r <= 40; r++){

	for(unsigned int s = 1; s <= 40; s++){

	map.first(r, s) = -2.0 + counter1;
    map.second(r, s) = -2.0 + counter2;
	counter1 = counter1+0.1;
                            }


counter1 = 0;
counter2 = counter2+0.1;

}

}



double fn (double r, double s, double rl, double sl) {

/* "Schreiben Sie eine Funktion, die für zwei gegebene Neuronen das Ergebnis der Nachbarschaftsfunktion
fN innerhalb der Kohonen-Karte berechnet" -> */

	double arg = -((r-rl)*(r-rl)+(s-sl)*(s-sl))/fnPara;

	return exp(arg);

}


void learnUpdate(double x1, double x2, double x3, double nu){


	double m = 0;
	unsigned int wr = 1;
	unsigned int ws = 1;
	double neighbor = 0;

	vector< vector <double>> min;

/*

"Schreiben Sie eine Funktion, die für einen gegebenen Eingabevektor x das Best Matching Neuron
z(r,s) bestimmt. Geben Sie dabei das Indexpaar (r;s) des Best Matching Neurons in der Kohonen-
Karte zurück" ->

 */
	
	
	for(unsigned int r = 1; r <= map.R; r++){

		vector<double> row;

		for(unsigned int s = 1; s <= map.S; s++){
			m = ((map.first(r, s)-x1)*(map.first(r, s)-x1)) + ((map.second(r, s)-x2)*(map.second(r, s)-x2)) + ((map.third(r, s)-x3)*(map.third(r, s)-x3));

			row.push_back(m);
			}


			min.push_back(row);

	}



	m = min[0][0];




	for(unsigned int r = 0; r < min.size(); r++){

		for(unsigned int s = 0; s < min[r].size(); s++){

			if(min[r][s] < m){
			m = min[r][s];

			wr = r+1; // <- "Geben Sie dabei das Indexpaar (r;s) des Best Matching Neurons zurück"
			ws = s+1; // <-

			}

		}

	} 

/* es ist jetzt zwar keine seperate Funktion vorhanden die das winner-neuron festlegt, aber auch gar nicht nötig und würde damit das Programm nur verlangsamen */


// "Implementieren Sie die Kohonen-Lernregel für einen gegebenen Eingabevektor x" ->


	for(unsigned int zr = 1; zr <= map.R; zr++){  
		for(unsigned int zs = 1; zs <= map.S; zs++){
	  
	  neighbor = fn((double)wr, (double)ws, (double)zr, (double)zs);// Abstand des Neurons vom winner-Neuron
// "Die Zentren der Neuronen werden durch die Kohonen-Lernregel adaptiert" ->

	  map.first(zr, zs) += nu * neighbor * (x1-map.first(zr, zs));
	  map.second(zr, zs) += nu * neighbor * (x2-map.second(zr, zs));
	  map.third(zr, zs) += nu * neighbor * (x3-map.third(zr, zs));

		                                  }

									  }


}


/* "Als Eingabeverteilung werden Vektoren x verwendet, deren ersten beiden Komponenten x1 und x2
zufällig aus dem Bereich [-2; 2] gezogen werden, während die dritte Komponente x3 wie folgt berechnet
wird" ->
*/

double f_X3(double x1, double x2) {


double x3 = 0;
double r = 2.0/3.0;
double tau;
double rho;
double distance;
double epsilon;

knn::init(); 
epsilon = knn::randomFromInterval(-0.1,0.1);



distance = (x1*x1)+(x2*x2);
distance = sqrt(distance);

if(r < distance && distance < 3*r){

tau = atan2(x2,x1);
rho = (x1/cos(tau))-(2.0*r);
rho = rho/r;
rho = acos(rho);

x3 = (r*sin(rho)) + epsilon;
}

else{
x3 = epsilon;}


return x3;

}

/*
"Führen Sie die Kohonen-Lernregel iterativ für viele mit der Funktion aus e)" ->[der Grund warum wir den x-Vektor nicht mit Datenstruktur und Funktion implementieren ist, das es vorher damit Probleme gab (immer die selbe Zahl als Ergebnis) und es eigentlich gar nicht notwendig ist da man zu einer Zeit nur einen Vektor braucht den man mit simplen Variablen direkt und neu erzeugen kann] "und erzeugten Eingabevek-
toren x durch. Geben Sie nach 1000, 10000, 100000 und 500000 Schritten die Zentrumsvektoren
μ(r,s) der Neuronen in eine Datei aus und plotten Sie diese." 

->

*/

void formMap(unsigned int iteration, double learnRate){

double x1 = 0;
double x2 = 0;
double x3 = 0;


//-----------------------------------------

{ //Block

ofstream initial ("0iterations.txt");

for (unsigned int r = 1; r <= map.R; ++r){
for (unsigned int s = 1; s <= map.S; ++s){
initial << map.first(r, s) << "\t" << map.second(r, s) << "\t" << map.third(r, s) << "\n";
}

}

initial.close();

}//Block end

//Start Process-------------------------------------------------------------

for(int i = 0; i < iteration; i++) {

knn::init();
x1 = knn::randomFromInterval(-2.0,2.0); // <- "zufällig aus dem Bereich [-2; 2]"
x2 = knn::randomFromInterval(-2.0,2.0);
x3 = f_X3(x1,x2);

learnUpdate(x1,x2,x3,learnRate);

cout << "Iteration: " << i << "\n";

if(i==999)
break;

}


/////////////////////////////////////////////////////////////////

{ //Block

ofstream tausend ("1.000iterations.txt"); // plotting the file

for (unsigned int r = 1; r <= map.R; ++r){
for (unsigned int s = 1; s <= map.S; ++s){
tausend << map.first(r, s) << "\t" << map.second(r, s) << "\t" << map.third(r, s) << "\n";
}

}

tausend.close();

}//Block end

//////////////////////////////////////////////////////////

for(int i = 1000; i < iteration; i++) {

knn::init();
x1 = knn::randomFromInterval(-2.0,2.0);
x2 = knn::randomFromInterval(-2.0,2.0);
x3 = f_X3(x1,x2);

learnUpdate(x1,x2,x3,learnRate);

cout << "Iteration: " << i << "\n";

if(i==9999)
break;

}

///////////////////////////////////////////////////////////////////

{ //Block

ofstream ten ("10.000iterations.txt");

for (unsigned int r = 1; r <= map.R; ++r){
for (unsigned int s = 1; s <= map.S; ++s){
ten << map.first(r, s) << "\t" << map.second(r, s) << "\t" << map.third(r, s) << "\n";
}

}

ten.close();

}

for(int i = 10000; i < iteration; i++) {

knn::init();
x1 = knn::randomFromInterval(-2.0,2.0);
x2 = knn::randomFromInterval(-2.0,2.0);
x3 = f_X3(x1,x2);

learnUpdate(x1,x2,x3,learnRate);

cout << "Iteration: " << i << "\n";

if(i==99999)
break;

}

///////////////////////////////////////////////////////////////////////////

{ //Block

ofstream hundert ("100.000iterations.txt");

for(unsigned int r = 1; r <= map.R; ++r){
for(unsigned int s = 1; s <= map.S; ++s){
hundert << map.first(r, s) << "\t" << map.second(r, s) << "\t" << map.third(r, s) << "\n";
}

}

hundert.close();

}

///////////////////////////////////////////////////////


for(int i = 100000; i < iteration; i++) {

knn::init();
x1 = knn::randomFromInterval(-2.0,2.0);
x2 = knn::randomFromInterval(-2.0,2.0);
x3 = f_X3(x1,x2);

learnUpdate(x1,x2,x3,learnRate);

cout << "Iteration: " << i << "\n";

/* if(i==499999)
break; */

}

//------------------------------------------------------------

ofstream endplot ("RESULT.txt");

for (unsigned int r = 1; r <= map.R; ++r){
for (unsigned int s = 1; s <= map.S; ++s){
endplot << map.first(r, s) << "\t" << map.second(r, s) << "\t" << map.third(r, s) << "\n";
}

}

endplot.close();

}



};


}
