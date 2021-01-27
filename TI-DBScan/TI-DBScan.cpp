// DAMI-DBSCAN.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono> 
#include <utility>
#include <algorithm>
//TODO: 1. make it work!2. correct border-core assignment 3.merge with DBSCAN
using namespace std;
int status;
const int NOT_CLASSIFIED = -1;
const int NOISE = -2;
class Point {
public:
    int index;
    double x, y;
    int noOfPoints;
    int clusterNo;
    int distances = 0;
    int type = 2;
    double distance = 0; //distance to ref
    double getDistance(const Point& point) {
        double distance = sqrt((x - point.x) * (x - point.x) + (y - point.y) * (y - point.y));
        //cout << "distance " << distance << endl;
        return distance;
    }
};
class TiDbScan {
public:
    int minPoints, clusterInx;
    double eps;
    vector<Point> points;
    vector<Point> points_sorted;
    int size;
    vector<vector<Point> > neighbours;
    vector<vector<int> > cluster;
    vector<int> noise;
    int ref;

    TiDbScan(double eps, int minPoints, vector<Point> points) {
        this->eps = eps;
        this->minPoints = minPoints;
        this->points = points;
        this->size = (int)points.size();
        neighbours.resize(size);
        noise.resize(size);
        this->points_sorted=points;//.resize(size);
        this->clusterInx = -1;
        //  this->corePoints.resize(size);
    }

    void run() {
        auto start = chrono::high_resolution_clock::now();
        //set reference point
        ref = points[size - 1].index - 1;
        //find distance of every point to the ref
        for (int i = 0; i < size; i++) {
            points[i].getDistance(points[ref]);
        }
        //sort points by distance to reference
        std::sort(points_sorted.begin(), points_sorted.end(), [](const Point& lhs, const Point& rhs)
            {
                return lhs.distance < rhs.distance;
            });
        //findNeighbours();
        auto end = chrono::high_resolution_clock::now();
        double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();

        time_taken *= 1e-9;
        //cout << time_taken;
        for (int i = 0; i < size; i++) {
            if (points[i].clusterNo != NOT_CLASSIFIED) continue;

            //check if expand cluster
            if (isCore(i)) {
                TI_formCluster(i, ++clusterInx);
                // corePoints.push_back(i);
            }
            /*if (isBorder(i)) {
                points[i].type = 0;}

            //noise*/
            else {
                points[i].clusterNo = NOISE;
                points[i].type = -1;
                noise.push_back(i);
            }
        }
        //expand the cluster 
        /*
        cluster.resize(clusterInx + 1);
        for (int i = 0; i < size; i++) {
            if (points[i].clusterNo != NOISE) {
                // cout << points[i].clusterNo <<"; " <<clusterInx<<endl;
                cluster[points[i].clusterNo].push_back(i);
            }
        }*/
        /* for (size_t  i = 0; i < cluster.size(); i++) {

             for (size_t  j = 0; j < cluster[i].size(); j++) {
                 cout << "cluster:" << i << ": "<< cluster[i][j] << endl;
                 }
             }
         for (size_t j = 0; j < noise.size(); j++) {
             cout << "cluster noise: " << noise[j] << ": " << noise.size()<<   endl;
         }*/
         //writeOut();
       //  printNb();
        }
    
    //output results to file
    void writeOut() {
        ofstream outputfile("OUT-TI");
        outputfile << "index," << "x," << "y," << "type," << "distances," << "cluster" << endl;
        for (size_t i = 0; i < cluster.size(); i++) {
            for (size_t j = 0; j < cluster[i].size(); j++) {
                int pId = cluster[i][j];
                outputfile << pId << "," << points[pId].x << "," << points[pId].y << "," << points[pId].type << "," << points[pId].distances << "," << i << endl;
            }
        }

    }
    void printN() {
        ofstream output("OUT-big");
        for (int i = 0; i < size; i++) {

            output << points[i].index << "," << points[i].clusterNo << "," << points[i].x << "," << points[i].y << endl;

        }
    }

    void TI_formCluster(int now, int clust) {
        vector<Point> neighbours = TI_findNeighbours(now);
       /* if (neighbours.size() < MinPts) {
            points[i].clusterNo = NOISE;
            return false;
        }
        }*/
          points[now].clusterNo = clust;
          if (!isCore(now))
          {
              if (isBorder(now)) { points[now].type = 0; }
              return;
          }

          for (auto& next : neighbours[now]) {
              if (points[next].clusterNo != NOT_CLASSIFIED && points[next].clusterNo != NOISE) continue;
              TI_formCluster(next, clust);
          }
    }
    //find neighbourhood
    vector<int> TI_findNeighbours(int e) {
        vector<Point> fwdNeighbourhood = TI_fwdNeighbourhood(e);
        vector<Point> neighbourhood = TI_bwNeighbourhood(e);
        neighbourhood.insert(neighbourhood.end(), fwdNeighbourhood.begin(), fwdNeighbourhood.end());
        vector<int>neighbourhood_int;
        neighbourhood_int.resize(neighbourhood.size());
        for (size_t i = 0; i < neighbourhood.size();i++) {
            neighbourhood_int[(int)i].push_back(neighbourhood[i].index);
        }
        return  neighbourhood_int;
    }
    vector<Point> TI_fwdNeighbourhood(int e) {
        std::vector<Point>fw;
        //find index of element in the sorted data
        auto it = find_if(points_sorted.begin(), points_sorted.end(), [&e](const Point& obj) {return obj.index == e; });
        int ind;
           if (it != points_sorted.end())
           {
                // found element. it is an iterator to the first matching element.
                // if you really need the index, you can also get it:
                 ind = std::distance(points_sorted.begin(), it);
           }
            double threshold = eps + points_sorted[ind].distance;
            for (int i = ind + 1; i < points_sorted.size(); i++) {
                // if q.dist > forwardThreshold then // q.dist � p.dist > Eps?
                // break;
                // endif
                if (points_sorted[i].distance > threshold) {
                    break;
                }
                // if Distance2(q, p) <= Eps then
                // append q to seeds;
                // endif
                if (points_sorted[ind].getDistance(points_sorted[i]) <= (eps)) {
                    fw.push_back(points_sorted[i]);
                }
            }
            return fw;
    }
    vector<Point> TI_bwNeighbourhood(int e) {
        vector<Point>bw;
        //find index of element in the sorted data
        auto it = find_if(points_sorted.begin(), points_sorted.end(), [&e](const Point& obj) {return obj.index == e; });
        int ind;
        if (it != points_sorted.end())
        {
            // found element. it is an iterator to the first matching element.
            // if you really need the index, you can also get it:
            ind = std::distance(points_sorted.begin(), it);
        }
        double threshold = eps - points_sorted[ind].distance;
        for (int i = ind - 1; i >=0; i--) {
            // if q.dist > forwardThreshold then 
            // break;
            // endif
            if (points_sorted[i].distance < threshold) {
                break;
            }
            // if Distance2(q, p) <= Eps then
            // append q to seeds;
            // endif
            if (points_sorted[ind].getDistance(points_sorted[i]) <= (eps)) {
                bw.push_back(points_sorted[i]);
            }
        }
        return bw;
    }
    // check if neighbourhood satisfies the minpoins requirement
    bool isCore(int index) {
        points[index].type = 1;
        return points[index].noOfPoints >= minPoints;
    }
    bool isBorder(int ind) {
        for (size_t i = 1; i < neighbours[ind].size(); i++) {
            if (isCore(neighbours[ind][i])) {
                return true;
            }
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
            fin >> index >> x >> y;
            /* while (getline(fin, element, ','))
             {

             }*/
            points.push_back({ index,x,y,0, NOT_CLASSIFIED });
        }
        points.pop_back();


        for (auto i = points.begin(); i != points.end(); i++)
        {
            cout << i->x << ' ' << i->y << ' ' << i->index << endl;

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
    std::string path = "file";
    InputReader input = InputReader(path);
    TiDbScan alg = TiDbScan(3, 1, input.getPoints());
    alg.run();
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
