#include "som.h"
#include<iostream>

/*

Michael Blachetta


"die Zentrumsvektoren
μ(r,s) der Neuronen in eine Datei aus und plotten Sie diese" ->
Im Ordner Outputs sind die Plots zu sehen, wir haben jeweils ein grid von fn mit 2.0 und eines mit 10.0 

zu Aufgabe f):

auf dem Gittergraphen der beiden Plots fn 2 und fn 10 ist deutlich zu sehen das der angedeutete Thorus
wenn bei fn mit 10 dividiert wird am äußeren Rand viel flacher/ebener ist als umgekehrt. Es fehlt eine gleichmäßige 
Wölbung vom Rand nach oben was die runde Form des Elipsoids ausmacht die bei fn 2, wenn auch mit großen Beulen, zu erkennen ist. 
Die innere Einbeulung in der Mitte des Thorus ist
etwas stärker ausgeprägt als beim anderen aber sonst gibt es Stellen die sehr stark von der Form abweichen u. a. am Rand.

der Faktor mit dem der euklidische Abstand im Exponentialteil dividiert wird ist ein Wert der das Sigma des Nachbarschafts und
Adaptionsradius in der Kohonenkarte vergrößert oder verkleinert. Je größer der Divisor, desto größer der Radius und desto mehr 
Neuronen werden in einem Lernzyklus von der Eingabe beeinflusst und verändert.

Bei fn 10 werden nun viel mehr Neuronen von einigen Eingaben hingezogen als nötig oder sinnvoll, 
da sie dann an anderer Stelle fehlen weil das Gitter nur eine begrenzte Anzahl von Neuronen hat.
Es sollten für jeden Input etwa gleichmäßig viele Punkte zur verfügung stehen 
die in die Stelle der Form integriert werden kann. 
In einigen Fällen jedoch kann der Lernprozess dadurch beschleunigt werden.

*/


using namespace std;

int main(void){

double n;
double sig;
unsigned int it;


cout << "SELF ORGANIZING MAP Thorus Simulator (by Michael Blachetta and Negah Moradi) \n\n\n";

cout << "enter learning interval as x.x \n"; // <- "wählen Sie als Lernrate nu = 0.1." sie können wählen was immer sie wünschen
cin >> n;

cout << "\nchoose map neighborhood distance factor as x.x \n"; // <- "Was ändert sich daran, wenn sie im Exponenten von fN durch 10 und nicht durch 2 teilen und warum", probieren sie auch mal fn = 5 aus
cin >> sig;

cout << "\nselect number of iterations (as integer) \n"; // <- "1000, 10000, 100000 und 500000 Schritte", sie können auch 1 Million wählen"
cin >> it;

cout << "\n\n initialising neuronal grid... \n\n";

som::neuroGrid thorus = som::neuroGrid(sig);

cout << "start map forming process... \n";

thorus.formMap(it,n); // <- "Führen sie die Kohonen-Lernregel iterativ für viele mit der Funktion aus e) erzeugten Eingabevektoren x durch"

//thorus.formMap(500000,0.1);

return 0;


}
