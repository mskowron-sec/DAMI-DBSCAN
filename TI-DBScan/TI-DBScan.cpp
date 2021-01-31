// DAMI-DBSCAN.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono> 
#include <utility>
#include <algorithm>
//TODO: 2. correct border-core assignment 3.merge with DBSCAN
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
    double distance_to_ref = 0;
    int distances = 0;
    int type = 2;
    double distance = 0; //distance to ref
    double getDistance( Point& point) {
        double distance = sqrt((x - point.x) * (x - point.x) + (y - point.y) * (y - point.y));
        this->distances++;
        point.distances++;
       // cout << "distance " << distance << endl;
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
    vector<vector<int> > neighbours;
    vector<vector<int> > cluster;
    Point ref;

    TiDbScan(double eps, int minPoints, vector<Point> points) {
        this->eps = eps;
        this->minPoints = minPoints;
        this->points = points;
        this->size = (int)points.size();
        neighbours.resize(size);
        this->points_sorted=points;//.resize(size);
        this->clusterInx = -1;
        Point ref = { 0,0,0,0, NOT_CLASSIFIED };
        //  this->corePoints.resize(size);
    }

    void run() {
        auto start = chrono::high_resolution_clock::now();
        
        //ref = points[size - 1].index - 1;
        //find distance of every point to the ref
        for (int i = 0; i < size; i++) {
            points_sorted[i].distance = points[i].getDistance(ref);
            points_sorted[i].distances++;
        }
        auto end = chrono::high_resolution_clock::now();
        double ref_time =
            chrono::duration_cast<chrono::microseconds>(end - start).count();
        //sort points by distance to reference
         start = chrono::high_resolution_clock::now();
        std::sort(points_sorted.begin(), points_sorted.end(), [](const Point& lhs, const Point& rhs)
            {
                return lhs.distance < rhs.distance;
            });
        //findNeighbours();
         end = chrono::high_resolution_clock::now();
        double sort_time =
            chrono::duration_cast<chrono::microseconds>(end - start).count();
        //find neighbourhood for every point in the dataset
        start = chrono::high_resolution_clock::now();
        for (int i = 0; i < size; i++) {
            neighbours[i] = TI_findNeighbours(i);
            points[i].noOfPoints = neighbours[i].size();
        }
        end = chrono::high_resolution_clock::now();
        double findneighbours_time =chrono::duration_cast<chrono::microseconds>(end - start).count();
        //cout << time_taken;
        start = chrono::high_resolution_clock::now();
        for (int i = 0; i < size; i++) {
            if (points[i].clusterNo != NOT_CLASSIFIED) {
                if (points[i].type != 1 && points[i].type != -2) { points[i].type = 0; }
                 continue;
            }
        
            //check if expand cluster
            if (isCore(i)) {
                points[i].type = 1;
                clusterInx++;
                TI_formCluster(i, clusterInx);
                // corePoints.push_back(i);
            }
            /*if (isBorder(i)) {
                }

            //noise*/
            else {
                points[i].clusterNo = NOISE;
                points[i].type = -1;
            }
        }
        end = chrono::high_resolution_clock::now();
        double classification_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
        start  = chrono::high_resolution_clock::now();
        //ref_time, sort_time,findneighbours_time, classification_time
        writeOutput();
        end = chrono::high_resolution_clock::now();
        double write_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout << write_time << endl;
        printNb();
        }
    //double reft, double sort, double neig, double clas
    //output results to file
    void writeOutput() {
        ofstream outputfti("TI-OUT");
       // ofstream statsti("TI-STAT");
        outputfti << "index," << "x," << "y," << "type," << "distances," << "cluster" << endl;
        for (int pId = 0; pId < size; pId++) {
            cout << points[pId].distance << endl;
            outputfti << points[pId].index << "," << points[pId].x << "," << points[pId].y << "," << points[pId].type << "," << points[pId].distances << "," << points[pId].clusterNo << endl;

        }
        //statsti << "eps=" << eps << endl << "minPoints=" << minPoints << endl << "pointsNo=" << size << endl << "clusters=" << clusterInx+1 << endl << "Time of distance calc: : " << reft << endl<< "Time of sorting: " << sort << endl<< "FindNeighbours time in microsec = " << neig << endl << "Time of classification: " << clas << endl;
    }
    void printN() {
        ofstream output("OUT-big");
        for (int i = 0; i < size; i++) {

            output << points[i].index << "," << points[i].clusterNo << "," << points[i].x << "," << points[i].y << endl;

        }
    }
    void printNb() {
        ofstream output("neighbours");
        for (int i = 0; i < size; i++) {
            output << i << endl;
            for (size_t j = 0; j < neighbours[i].size(); j++) {
                //int pId = cluster[i][j];
                output << neighbours[i][j] << ",";// << points[pId].y << "," << points[pId].type << "," << points[pId].distances << "," << i << endl;
            }
            output << endl;
        }
    }
    void TI_formCluster(int now, int c) {
        points[now].clusterNo = c;
        if (!isCore(now))
        {
           if (points[now].type != 1 && points[now].type != -2) { points[now].type = 0; }
            return;
        }
        else { points[now].type = 1; }

        for (int i = 0; i < (int)neighbours[now].size(); i++) {
            int next = neighbours[now][i];
            if (points[next].clusterNo != NOT_CLASSIFIED && points[next].clusterNo != NOISE) {   continue; }
            TI_formCluster(next, c);
        }
    }
    //find neighbourhood
    vector<int> TI_findNeighbours(int e) {
        int ind = 0;
        //find element in the sorted table
        auto it = find_if(points_sorted.begin(), points_sorted.end(), [&e](const Point& obj) {return obj.index == e; });

        if (it != points_sorted.end())
        {
            // found element. it is an iterator to the first matching element.
            // if you really need the index, you can also get it:
            ind = std::distance(points_sorted.begin(), it);
        }
        cout << ind << endl;
        vector<int> fwdNeighbourhood = TI_fwdNeighbourhood(ind,e);
        vector<int> neighbourhood = TI_bwNeighbourhood(ind,e);
        neighbourhood.insert(neighbourhood.end(), fwdNeighbourhood.begin(), fwdNeighbourhood.end());
       /* vector<int>neighbourhood_int;
        neighbourhood_int.resize(neighbourhood.size());
        for (size_t i = 0; i < neighbourhood.size();i++) {
            neighbourhood_int[(int)i].push_back(neighbourhood[i].index);
        }*/
        return  neighbourhood;
    }
    vector<int> TI_fwdNeighbourhood(int e,int old) {
        std::vector<int>fw;
            double threshold = eps + points_sorted[e].distance;
            for (int i = e + 1; i < points_sorted.size(); i++) {
                // if q.dist > forwardThreshold then // q.dist � p.dist > Eps?
                // break;
                // endif
                if (points_sorted[i].distance > threshold) {
                    break;
                }
                // if Distance2(q, p) <= Eps then
                // append q to seeds;
                // endif
                double dist = points_sorted[e].getDistance(points_sorted[i]);
                points_sorted[i].distances++;
                points_sorted[e].distances++;
                if (dist <= (eps)) {
                    fw.push_back(points_sorted[i].index);
                }
            }
            return fw;
    }
    vector<int> TI_bwNeighbourhood(int e,int old) {
        vector<int>bw;
        //find index of element in the sorted data
        auto it = find_if(points_sorted.begin(), points_sorted.end(), [&e](const Point& obj) {return obj.index == e; });
        int ind;
        if (it != points_sorted.end())
        {
            // found element. it is an iterator to the first matching element.
            // if you really need the index, you can also get it:
            ind = std::distance(points_sorted.begin(), it);
        }
        double threshold = eps - points_sorted[e].distance;
        for (int i = e - 1; i >=0; i--) {
            // if q.dist > forwardThreshold then 
            // break;
            // endif
            if (points_sorted[i].distance < threshold) {
                break;
            }
            // if Distance2(q, p) <= Eps then 
            // append q to seeds;
            // endif
            points_sorted[i].distances++;
            points[old].distances++;
            if (points_sorted[e].getDistance(points_sorted[i]) <= (eps)) {
                bw.push_back(points_sorted[i].index);
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
        for (size_t i = 1; i < neighbours[(int)ind].size(); i++) {
            if (isCore(neighbours[(int)ind][(int)i])) {
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

int main(int argc, char* argv[])
{
    
    std::string path = argv[1];
    InputReader input = InputReader(path);

    double epsilon = atoi(argv[2]);
     int minp = atoi(argv[3]);
    TiDbScan alg = TiDbScan(epsilon, minp, input.getPoints());
   /*
     std::string path = "file"; 
    InputReader input = InputReader(path);
    
    TiDbScan alg = TiDbScan(2, 3, input.getPoints());
    alg.run();
    */
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
