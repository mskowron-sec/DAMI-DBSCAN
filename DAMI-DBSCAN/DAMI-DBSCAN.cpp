// DAMI-DBSCAN.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono> 
//TODO: cmd arguments, more dimensions for points, calc number and core-border corrections 
using namespace std;
//int status;
const int NOT_CLASSIFIED = -1;
const int NOISE = -2;
class Point {
public:
    int index;
    double x, y;
    int noOfPoints;
    int clusterNo;
    int distances =0;
    int type = 0;
    double getDistance(const Point& point) {
        double distance=  sqrt((x - point.x) * (x - point.x) + (y - point.y) * (y - point.y));
        //cout << "distance " << distance << endl;
        return distance;
    }
};
class DbScan {
public:
    int minPoints, clusterInx;
    double eps;
    vector<Point> points;
    int size;
    vector<vector<int> > neighbours;
    vector<vector<int> > cluster;
    vector<int> noise;

    DbScan(double eps, int minPoints, vector<Point> points) {    
        this->eps = eps;
        this->minPoints = minPoints;
        this->points = points;
        this->size = (int)points.size();
        neighbours.resize(size);
        noise.resize(size);
        this->clusterInx = -1;
      //  this->corePoints.resize(size);
    }
    void run() {
        auto nstart = chrono::high_resolution_clock::now();
        findNeighbours();
        auto nend = chrono::high_resolution_clock::now();
        auto time_neighb = chrono::duration_cast<chrono::microseconds>(nend - nstart).count();;
       for (int i = 0; i < size; i++) {
            if (points[i].clusterNo != NOT_CLASSIFIED) continue;
            if (isCore(i)) {
                points[i].type = 1;
                formCluster(i, ++clusterInx);
               // corePoints.push_back(i);
            }
          //  if (isBorder(i)) {
            //    points[i].type = 0;}
            
            //noise
            else {
                points[i].clusterNo = NOISE;
                points[i].type = -1;
                noise.push_back(i);
            }
        }
       auto aend= chrono::high_resolution_clock::now();
       auto time_all = chrono::duration_cast<chrono::microseconds>(aend - nstart).count();
        //cluster structures
        /*cluster.resize(clusterInx + 1);
        for (int i = 0; i < size; i++) {
            if (points[i].clusterNo != NOISE) {
               // cout << points[i].clusterNo <<"; " <<clusterInx<<endl;
                cluster[points[i].clusterNo].push_back(i);
            }
        }*/

        writeOutput(time_neighb, time_all);
 
        printNb();

    }
    //output results to file
    void writeOutClusters() {
        ofstream outputfile("OUT2");
       // ofstream stats("STAT");
        outputfile <<"index,"<<"x,"<<"y,"<<"type,"<<"distances,"<<"cluster"  <<endl ;
        for (size_t i = 0; i < cluster.size(); i++) {
            for (size_t j = 0; j < cluster[i].size(); j++) {
                int pId = cluster[i][j];
                outputfile << pId << "," <<points[pId].x <<","<< points[pId].y << ","<< points[pId].type << ","<< points[pId].distances <<"," << i << endl;
            }
        }
    }
    void writeOutput(long time_neighb, long  time_all) {
        ofstream outputf("OUT");
        ofstream stats("STAT");
        outputf << "index," << "x," << "y," << "type," << "distances," << "cluster" << endl;
        for (int pId = 0; pId < size; pId++) {
            
                outputf << points[pId].index << "," << points[pId].x << "," << points[pId].y << "," << points[pId].type << "," << points[pId].distances << "," << points[pId].clusterNo << endl;
            
        }
        stats << "eps=" << eps << endl << "minPoints=" << minPoints << endl << "pointsNo=" << size << endl << "clusters=" << clusterInx+1 << endl << "FindNeighbours time in microsec = " << time_neighb << endl << "Time Overall in microsec= " << time_all << endl;;
    }
    void printN() {
        ofstream output("OUT-big");
        for (int i = 0; i < size; i++) {

                output << points[i].index << "," << points[i].clusterNo <<  "," << points[i].x << "," << points[i].y << endl;
            
        }
    }
    void printNb() {
        ofstream output("neighbours");
        for (int i = 0; i < size; i++) {
            output << i << endl;
            for (size_t j = 0; j < neighbours[i].size(); j++) {
                //int pId = cluster[i][j];
                output  << neighbours[i][j] << ",";// << points[pId].y << "," << points[pId].type << "," << points[pId].distances << "," << i << endl;
            }
            output << endl;
        }
    }
    //find neighbourhood
    void findNeighbours() {
        
        for (int i = 0; i < size; i++) {
           // neighbours[i].push_back(i);
           // points[i].noOfPoints++;
            for (int j = 0; j < size; j++) {
                if (i == j) continue;
                //if within epsilon radius
                double dist = points[i].getDistance(points[j]);
                points[j].distances++;
                points[i].distances++;
                if (dist <= eps) {
                    points[i].noOfPoints++;
                    //push to neighbourhood
                    neighbours[i].push_back(j);
                }
            }
        }
    }
    // check if neighbourhood satisfies the minpoins requirement
    bool isCore(int index) {
        
        return points[index].noOfPoints >= minPoints;
    }
    bool isBorder(int ind) {
        for (size_t i = 1; i < neighbours[ind].size(); i++) {
            if (isCore(neighbours[ind][i])) {
                return true;
            }
        }
    }
    void formCluster(int now, int c) {
        points[now].clusterNo = c;
        if (!isCore(now)) 
        {
            if (isBorder(now)) { points[now].type = 0; }
            return;
        }
        else{ points[now].type = 1; }
        
        for (int i=0;i<(int)neighbours[now].size();i++) {
            int next = neighbours[now][i];
            if (points[next].clusterNo != NOT_CLASSIFIED && points[next].clusterNo != NOISE) {
                points[now].type = 0;  continue;
             }
            formCluster(next, c);
        }
    }
    vector<vector<int> > getCluster() {
        return cluster;
    }
};
class InputReader {
private:
    ifstream fin;
    vector<Point> points;
public:
    InputReader(string filename) {
        fin.open(filename.c_str());
        if (!fin) {
              cout << filename << " could not be read - file not found!\n";
             exit(0);
        }
        parse(); 
    }
    void parse() {
        int index;
        double x, y;
        string element;
        while (!fin.eof()) {
            fin >>index >> x >> y;
           /* while (getline(fin, element, ','))
            {
               
            }*/
            points.push_back({ index,x,y,0, NOT_CLASSIFIED});
        }
        points.pop_back();
         
         
        for (auto i = points.begin(); i != points.end(); i++)
        {
            cout << i->x<<' '<< i->y <<' ' <<i->index << endl;

        }
    }
    vector<Point> getPoints() {
        return points;
    }
};
//

int main(int argc, char* argv[])
{
    std::string path = argv[1];
    InputReader input = InputReader(path);
    
    double epsilon = atoi(argv[2]);
   // cout << epsilon << endl;

    int minp = atoi(argv[3]);
    //cout <<"minpts:"<< minp << endl;
    DbScan alg = DbScan(epsilon, minp, input.getPoints());
    alg.run();
    /*std::string path = "file";
    InputReader input = InputReader(path);

    DbScan alg = DbScan(2, 3, input.getPoints());
    alg.run();*/
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
