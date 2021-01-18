// DAMI-DBSCAN.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
using namespace std;
int status;
const int NOT_CLASSIFIED = -1;
class Point {
public:
    double x, y;
    int noOfPoints;
    int clusterNo;
    int label;
    double getDistance(const Point& point) {
        return sqrt((x - point.x) * (x - point.x) + (y - point.y) * (y - point.y));
    }
};
class InputReader {
private:
    ifstream fin;
    vector<Point> points;
    enum Status { NOT_VISITED = -1, BORDER, NOISE, CORE };
public:
    InputReader(string filename) {
        fin.open(filename);
        if (!fin) {
            cout << filename << "Error!\n";
            exit(0);
        }
        parse();
    }
    void parse() {
        int label;
        double x, y;
        string element;
        while (!fin.eof()) {
            //fin >> x>>" ">>y>>"," >> label;
            while (getline(fin, element, ','))
            {
               
            }
            points.push_back({ x,y,0, NOT_CLASSIFIED, label});
        }
        points.pop_back();
         
         
        for (auto i = points.begin(); i != points.end(); i++)
        {
            cout << i->x<<' '<< i->y <<' ' <<i->label << endl;

        }
    }
    vector<Point> getPoints() {
        return points;
    }
};
//Point class - 2D
//TODO N-Dimensions

int main()
{
    std::cout << "Hello World!\n";
    std::string path = "test.arff";
    InputReader input = InputReader(path);
   // std::cout << input.getPoints()
    
}

// Uruchomienie programu: Ctrl + F5 lub menu Debugowanie > Uruchom bez debugowania
// Debugowanie programu: F5 lub menu Debugowanie > Rozpocznij debugowanie

// Porady dotyczące rozpoczynania pracy:
//   1. Użyj okna Eksploratora rozwiązań, aby dodać pliki i zarządzać nimi
//   2. Użyj okna programu Team Explorer, aby nawiązać połączenie z kontrolą źródła
//   3. Użyj okna Dane wyjściowe, aby sprawdzić dane wyjściowe kompilacji i inne komunikaty
//   4. Użyj okna Lista błędów, aby zobaczyć błędy
//   5. Wybierz pozycję Projekt > Dodaj nowy element, aby utworzyć nowe pliki kodu, lub wybierz pozycję Projekt > Dodaj istniejący element, aby dodać istniejące pliku kodu do projektu
//   6. Aby w przyszłości ponownie otworzyć ten projekt, przejdź do pozycji Plik > Otwórz > Projekt i wybierz plik sln
