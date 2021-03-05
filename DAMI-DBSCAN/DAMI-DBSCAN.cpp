// DAMI-DBSCAN.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono> 
#include <sstream>
#include <cmath>
#include <assert.h>

using namespace std;
//int status;
const int NOT_CLASSIFIED = -2;
const int NOISE = -1;
class Point {
public:
    int index;
    vector<double>coordinates;
    // double x, y;
    int noOfNeighbours;
    int clusterNo;
    int label;
    int noDistances = 0;
    int type = 0;
    int visited = 0;
    double getEuclDistance2D(Point& point, int dim) {
        point.noDistances++;
        this->noDistances++;
        double distance = sqrt(pow(this->coordinates[0] - point.coordinates[0], 2) + pow(this->coordinates[1] - point.coordinates[1], 2));
        //cout << "distance " << distance << endl;
        return distance;
    }
    /// <summary>
    /// function to calculate Minkowski distance for N-dimensional vector, l-order
    /// </summary>
    /// <param name="point"></param>
    /// <param name="dim"></param>
    /// <param name="l"></param>
    /// <returns>distance</returns>
    double getDistanceN(Point& point, int dim, double l = 2.0) {

        double dim_sum = 0;
        for (int i = 0; i < dim; i++) {
            dim_sum += pow(coordinates[i] - point.coordinates[i], l);
        }
        double power = 1.0 / l;
        // cout << power << endl;
        double distance = pow(dim_sum, power);
        // cout << "distance " << distance << endl;
        return distance;
    };
};
class DbScan {
public:
    int minPoints, clusterInx;
    double eps;
    vector<Point> points;
    int size;
    vector<vector<int> > neighbours;
    int dimensions;


    DbScan(double eps, int minPoints, vector<Point> points) {
        this->eps = eps;
        this->minPoints = minPoints;
        this->points = points;
        this->size = (int)points.size();
        neighbours.resize(size);
        this->clusterInx = -1;
        this->dimensions = points[0].coordinates.size();
    }
    void run() {
        auto nstart = chrono::high_resolution_clock::now();
        //find neighbours functions
        findNeighbours();
        auto nend = chrono::high_resolution_clock::now();
        auto time_neighb = chrono::duration_cast<chrono::microseconds>(nend - nstart).count();
        cout << "Find neighbours time: " << time_neighb << endl;
        //printNb();
        auto cl = chrono::high_resolution_clock::now();
        //clustering
        for (int i = 0; i < size; i++) {

            if (points[i].clusterNo != NOT_CLASSIFIED) continue;

            if (isCore(i)) {
                points[i].type = 1;
                ++clusterInx;
                points[i].clusterNo = clusterInx;
                points[i].visited = 1;
                for (int n = 0; n < neighbours[i].size(); n++) {
                    points[neighbours[i][n]].visited = 1;
                }
                //search among neighbors
                for (int j = 0; j < neighbours[i].size(); j++) {
                    int cur_n = neighbours[i][j];

                    if (points[cur_n].clusterNo == NOISE) { points[cur_n].clusterNo = clusterInx; points[cur_n].visited = 1; points[cur_n].type = 0; }
                    if (points[cur_n].clusterNo != NOT_CLASSIFIED) { continue; }
                    points[cur_n].clusterNo = clusterInx;
                    points[cur_n].visited = 1;
                    if (isCore(cur_n)) {
                        points[cur_n].type = 1;
                        for (int k = 0; k < neighbours[cur_n].size(); k++)
                        {
                            int cand_idx = neighbours[cur_n][k];
                            if (points[cand_idx].clusterNo == NOISE) { points[cand_idx].clusterNo = clusterInx; }
                            if (points[cand_idx].clusterNo == NOT_CLASSIFIED) {
                                if (points[cand_idx].visited == 0)
                                {
                                    neighbours[i].push_back(neighbours[cur_n][k]);
                                    points[cand_idx].visited = 1;
                                }
                            }
                        }
                    }
                }
            }
            //noise
            else {
                points[i].clusterNo = NOISE;
                points[i].type = -1;
            }
        }
        auto aend = chrono::high_resolution_clock::now();
        auto time_class = chrono::duration_cast<chrono::microseconds>(aend - cl).count();
        cout << "Clustering time: " << time_class << endl;
        auto time_all = chrono::duration_cast<chrono::microseconds>(aend - nstart).count();

        auto wstart = chrono::high_resolution_clock::now();
        writeOutput(time_neighb, time_class, time_all);
        auto all = chrono::high_resolution_clock::now();
        auto time_write = chrono::duration_cast<chrono::microseconds>(all - wstart).count();
        auto time_all_w = chrono::duration_cast<chrono::microseconds>(all - nstart).count();
        cout << "Time of  writing: " << time_write << " microsecond" << endl;
        cout << "Time of run() with writing: " << time_all_w << " microsecond" << endl;

        //printNb();

    }
    //output results to file
    void writeOutput(long time_neighb, long clas, long  time_all) {
        //variable to store sum of distance calculations
        double dist_sum = 0;
        //average of distance calc per point
        double dist_avg = 0;
        int rand_idx = 0;
        ofstream outputf("OUT");
        ofstream stats("STAT");
        outputf << "index," << "coordinates," << "type," << "distances," << "cluster" << endl;
        for (int pId = 0; pId < size; pId++) {
            dist_sum += points[pId].noDistances;
            if (points[pId].clusterNo == points[pId].label)
            {
                ++rand_idx;
            }
            outputf << points[pId].index << ",";
            for (int j = 0; j < dimensions; j++) {
                outputf << points[pId].coordinates[j] << ",";
            }
            outputf << points[pId].type << "," << points[pId].noDistances << "," << points[pId].clusterNo << endl;

        }
        double rand_index = rand_idx / (double)size;
        dist_avg = dist_sum / size;
        stats << "eps=" << eps << endl << "minPoints=" << minPoints << endl << "pointsNo=" << size << endl << "dimensions=" << dimensions << endl << "clusters=" << clusterInx + 1 << endl << "FindNeighbours time in microsec = " << time_neighb << endl << "clustering time= " << clas << endl << "Time Overall of run() in microsec= " << time_all << endl << "Distance calc per point=" << dist_avg << endl << "Rand index=" << rand_index << endl;
        cout << "eps=" << eps << ' ' << "minPts=" << minPoints << ' ' << "clusters=" << clusterInx + 1 << endl;
    }
    //debug function to verify correct neighbour assignment
    void printNb() {
        ofstream output("neighbours");
        for (int i = 0; i < size; i++) {
            output << i << endl;
            for (size_t j = 0; j < neighbours[i].size(); j++) {

                output << neighbours[i][j] << ",";
            }
            output << endl;
        }
    }
    //find neighbourhood for all points
    void findNeighbours() {

        for (int i = 0; i < size; i++) {

            for (int j = i + 1; j < size; j++) {
                //if within epsilon radius
                double dist = points[i].getDistanceN(points[j], dimensions);

                //points[j].noDistances++;
                points[i].noDistances++;
                if (dist <= eps) {
                    points[i].noOfNeighbours++;
                    //push to neighbourhood
                    neighbours[i].push_back(j);

                    // increment neighbour count
                    points[j].noOfNeighbours++;
                    neighbours[j].push_back(i);
                }
            }

            // add self
            points[i].noOfNeighbours++;
            assert(points[i].noOfNeighbours == neighbours[i].size() + 1);

            points[i].noOfNeighbours = neighbours[i].size() + 1;


        }
    }
    ///  check if neighbourhood satisfies the minpoins requirement
    bool isCore(int index) {
        if (points[index].noOfNeighbours >= minPoints) {
            return true;
        }
        else { return false; }
    }
    /// <summary>
    /// check if is a border point
    /// </summary>
    /// <param name="filename"></param>
    bool isBorder(int ind) {
        for (size_t i = 1; i < neighbours[ind].size(); i++) {
            if (isCore(neighbours[ind][i])) {
                return true;
            }
        }
    }
    //old recursive function to form clsters
    void formCluster(int now, int c) {

        points[now].clusterNo = c;
        if (!isCore(now))
        {
            if (isBorder(now)) { points[now].type = 0; }
            return;
        }
        else { points[now].type = 1; }

        for (int i = 0; i < (int)neighbours[now].size(); i++) {
            int next = neighbours[now][i];
            //
            if (points[next].clusterNo != NOT_CLASSIFIED) {
                continue;
            }
            formCluster(next, c);
        }
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
        int index = 0;
        int label;
        double x, y;
        string coordinates;

        while (!fin.eof()) {
            fin >> label >> coordinates;
            ///vector to store coordinates
            vector<double> elements;
            std::stringstream ss(coordinates);
            while (ss.good())
            {

                string substr;
                getline(ss, substr, ',');

                elements.push_back(atof(substr.c_str()));
            }
            points.push_back({ index,elements, 0, NOT_CLASSIFIED,label });
            index++;
        }

    }
    /// <summary>
    /// getter for input point list
    /// </summary>
    /// <param name="argc"></param>
    /// <param name="argv"></param>
    /// <returns> vector<Point >points</returns>
    vector<Point> getPoints() {
        return points;
    }
};
//

int main(int argc, char* argv[])
{

    auto tstart = chrono::high_resolution_clock::now();
    // std::string path = "simple";
    std::string path = argv[1];
    InputReader input = InputReader(path);
    double epsilon = atof(argv[2]);
    int minp = atoi(argv[3]);
    auto wend = chrono::high_resolution_clock::now();
    auto read_time = chrono::duration_cast<chrono::microseconds>(wend - tstart).count();
    cout << "Read time: " << read_time << " microseconds" << endl;
    DbScan alg = DbScan(epsilon, minp, input.getPoints());
    auto rstart = chrono::high_resolution_clock::now();
    alg.run();
    auto rend = chrono::high_resolution_clock::now();
    auto run_time = chrono::duration_cast<chrono::microseconds>(rend - rstart).count();
    auto stop = chrono::high_resolution_clock::now();
    auto time_whole = chrono::duration_cast<chrono::microseconds>(stop - tstart).count();
    cout << "Run fuction time (in main): " << run_time << " microseconds" << endl;

    cout << "Total Execution time: " << time_whole;
}

// Uruchomienie programu: Ctrl + F5 lub menu Debugowanie > Uruchom bez debugowania
// Debugowanie programu: F5 lub menu Debugowanie > Rozpocznij debugowanie
