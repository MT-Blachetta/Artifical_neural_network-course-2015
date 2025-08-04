# Artificial Neural Network Course 2015

Dieses Repository enthält drei voneinander unabhängige C++-Programme, die im Rahmen einer Lehrveranstaltung zu künstlichen neuronalen Netzen entstanden sind. Jedes Programm demonstriert einen anderen Aspekt der Thematik und kann separat kompiliert und ausgeführt werden.

## Verzeichnisstruktur
- `Polynomregression` – Approximiert Funktionen mit einer Polynomregression.
- `neuronal networks` – Interaktiver Simulator für ein mehrschichtiges Feedforward-Netz mit Backpropagation.
- `self-organizing map` – Umsetzung einer zweidimensionalen Kohonen-Karte zur Abbildung eines 3D-Torus.

## Theorie und Funktionsumfang
### Polynomregression
Das Programm erzeugt Trainings- und Testdaten für Sinus- bzw. Sinc-Funktionen und passt Polynome unterschiedlichen Grades mithilfe der kleinsten Quadrate an. Für jeden Grad werden der mittlere quadratische Fehler auf Trainings- und Testmenge berechnet.

### Feedforward-Netz
Der Simulator erstellt ein mehrschichtiges Netz mit sigmoiden Aktivierungen und trainiert es per Error-Backpropagation. Eingabe- und Zielvektoren werden automatisch generiert; Gewichte lassen sich zufällig initialisieren und über mehrere Epochen anpassen.

### Self‑Organizing Map
Die SOM besteht aus einem 40×40-Gitter, dessen Neuronen jeweils einen 3D-Zentrumsvektor besitzen. Eine Gaußsche Nachbarschaftsfunktion passt die Zentren iterativ an zufällige Eingabevektoren an und bildet so einen Torus im Eingaberaum.

## Benutzung
### Kompilieren
```bash
g++ -std=c++11 Polynomregression/main.cpp -o polyreg
g++ -std=c++11 'neuronal networks/main.cpp' -o network
g++ -std=c++11 'self-organizing map/main.cpp' -o som
```

### Polynomregression ausführen
```bash
./polyreg [P] [MMax]
```
- `P` – Anzahl der Trainingspunkte (Standard 11).
- `MMax` – maximaler Polynomgrad (Standard 19).

Ausgabe: `SinError.txt` und `SincError.txt` mit Zeilen `M Fehler_train Fehler_test`.

### Feedforward-Netz ausführen
```bash
./network
```
Das Programm fragt interaktiv nach:
1. Dateinamenssuffix für Ergebnisdateien.
2. Anzahl Eingabeneuronen.
3. Aufbau versteckter Layer (wiederholte Abfrage).
4. Größe der Trainings/Testmenge.
5. Grenzen zur zufälligen Gewichtsinitialisierung.
6. Lernrate und Anzahl Epochen.

Ausgabe: Fehlerverlauf in `meanErrors_<suffix>` sowie Funktionswerte in weiteren Dateien.

### Self‑Organizing Map ausführen
```bash
./som
```
Interaktive Eingaben:
1. Lernrate `ν`.
2. Faktor der Nachbarschaftsfunktion `σ` (teilt den euklidischen Abstand im Exponenten).
3. Anzahl Iterationen.

Ausgabe: Dateien wie `0iterations.txt`, `1.000iterations.txt`, … mit den Koordinaten der Neuronen nach bestimmten Iterationszahlen.

## Interpretation der Ausgaben
- In den Fehlerdateien repräsentiert jede Zeile den Fehler einer bestimmten Modellkomplexität; ein kleinerer Wert bedeutet bessere Approximation.
- Die erzeugten Koordinatendateien der SOM können mit Plotprogrammen visualisiert werden; sie zeigen die Anpassung des Gitters an die Torusoberfläche.

## Voraussetzungen
Benötigt wird ein C++11‑kompatibler Compiler. Für die optionale Visualisierung des Feedforward-Netzes wird `gnuplot-iostream` verwendet, welches eine funktionierende Gnuplot‑ und Boost‑Installation voraussetzt.
