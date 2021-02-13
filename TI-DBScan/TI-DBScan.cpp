// DAMI-DBSCAN.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono> 
#include <utility>
#include <algorithm>
#include <sstream>


using namespace std;
int status;
//flags for points not visited and noise point
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
    double distance = 0; //distance to ref

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
    }
};
class TiDbScan {
public:
    int minPoints, clusterInx;
    double eps;

    vector<Point> points_sorted;
    int size;
    vector<vector<int> > neighbours;
    Point ref;
    int dimensions;


    TiDbScan(double eps, int minPoints, vector<Point> points) {
        this->eps = eps;
        this->minPoints = minPoints;
        //       this->points = points;
        this->size = (int)points.size();
        neighbours.resize(size);
        this->points_sorted = points;//.resize(size);
        this->clusterInx = -1;
        this->dimensions = points[0].coordinates.size();


    }
    //"main"  function of the algorithm - finds reference distance, sorts points by this value, finds neighbour and finally cluster points
    void run() {
        auto start = chrono::high_resolution_clock::now();
        vector<double> zeros(dimensions, 0.0);
        Point ref = { 0,zeros,0,0, NOT_CLASSIFIED };

        //find distance of every point to the ref
        for (int i = 0; i < size; i++) {
            points_sorted[i].distance = points_sorted[i].getDistanceN(ref, dimensions);
        }
        auto end = chrono::high_resolution_clock::now();
        double ref_time =
            chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout << "Reference calc time:  " << ref_time << endl;
        //sort points by distance to reference
        start = chrono::high_resolution_clock::now();
        std::sort(points_sorted.begin(), points_sorted.end(), [](const Point& lhs, const Point& rhs)
            {
                return lhs.distance < rhs.distance;
            });

        end = chrono::high_resolution_clock::now();
        double sort_time =
            chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout << "Sorting time:  " << sort_time << endl;

        //find neighbourhood for every point in the dataset
        start = chrono::high_resolution_clock::now();
        for (int i = 0; i < size; i++) {
            neighbours[i] = TI_findNeighbours(i);
            points_sorted[i].noOfNeighbours = neighbours[i].size() + 1;
        }
        end = chrono::high_resolution_clock::now();
        double findneighbours_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout << "Find neighbours time: " << findneighbours_time << endl;
        auto cstart = chrono::high_resolution_clock::now();
        for (int i = 0; i < size; i++) {
            if (points_sorted[i].clusterNo != NOT_CLASSIFIED) {

                continue;
            }
            //check if new cluster should be created
            if (isCore(i)) {
                points_sorted[i].type = 1;
                clusterInx++;
                TI_formCluster(i, clusterInx);
            }
            //noise*/
            else {
                points_sorted[i].clusterNo = NOISE;
                points_sorted[i].type = -1;
            }
        }
        auto cend = chrono::high_resolution_clock::now();
        double classification_time = chrono::duration_cast<chrono::microseconds>(cend - cstart).count();
        cout << "Clustering time: " << classification_time << endl;
        start = chrono::high_resolution_clock::now();
        //
        writeOutput(ref_time, sort_time, findneighbours_time, classification_time);
        end = chrono::high_resolution_clock::now();
        double write_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout << "Write to file time=" << write_time << endl;
        // printNb();
    }
    //
    //output results to file
    void writeOutput(double reft, double sort, double neig, double clas) {
        //variable to store sum of distance calculations
        double dist_sum = 0;
        //average of distance calc per point
        double dist_avg = 0;
        int rand_idx = 0;
        ofstream outputfti("TI-OUT");
        ofstream statsti("TI-STAT");
        outputfti << "sorted index," << "index," << "coordinates[" << dimensions << "],type," << "distances," << "cluster" << endl;
        for (int i = 0; i < size; i++) {
            dist_sum += points_sorted[i].noDistances;
            if (points_sorted[i].clusterNo == points_sorted[i].label)
            {
                ++rand_idx;
            }
            outputfti << i << "," << points_sorted[i].index << ",";
            for (int j = 0; j < dimensions; j++) {
                outputfti << points_sorted[i].coordinates[j] << ",";
            }
            outputfti << points_sorted[i].type << "," << points_sorted[i].noDistances << "," << points_sorted[i].clusterNo << endl;
        }
        double rand_index = rand_idx / (double)size;
        dist_avg = dist_sum / size;
        statsti << "eps=" << eps << endl << "minPoints=" << minPoints << endl << "pointsNo=" << size << endl << "clusters=" << clusterInx + 1 << endl << "dimensions=" << dimensions << endl << "Time of distance calc: : " << reft << endl << "Time of sorting: " << sort << endl << "FindNeighbours time in microsec = " << neig << endl << "Time of clustering: " << clas << endl << "Avg calculations per point=" << dist_avg << endl << "Rand index=" << rand_index << endl;
        cout << "eps=" << eps << ' ' << "minPts=" << minPoints << ' ' << "clusters=" << clusterInx + 1 << endl;
    }

    void printNb() {
        ofstream output("neighbours");
        for (int i = 0; i < size; i++) {
            output << i << " " << points_sorted[i].index << endl;
            for (size_t j = 0; j < neighbours[i].size(); j++) {
                //print neighbours - real index in brackets
                output << neighbours[i][j] << "(" << points_sorted[neighbours[i][j]].index << ")" << ",";
            }
            output << endl;
        }
    }
    // function assigning point to clusters - will be updated to iterative in the next revision
    void TI_formCluster(int now, int c) {
        points_sorted[now].clusterNo = c;
        if (!isCore(now))
        {
            if (isBorder(now)) { points_sorted[now].type = 0; }
            return;
        }
        else { points_sorted[now].type = 1; }

        for (int i = 0; i < (int)neighbours[now].size(); i++) {
            int next = neighbours[now][i];
            if (points_sorted[next].clusterNo != NOT_CLASSIFIED && points_sorted[next].clusterNo != NOISE) { continue; }
            TI_formCluster(next, c);
        }
    }
    //find neighbourhood of a point with index e
    vector<int> TI_findNeighbours(int e) {
        vector<int> fwdNeighbourhood = TI_fwdNeighbourhood(e);
        vector<int> neighbourhood = TI_bwNeighbourhood(e);
        neighbourhood.insert(neighbourhood.end(), fwdNeighbourhood.begin(), fwdNeighbourhood.end());

        return  neighbourhood;
    }
    vector<int> TI_fwdNeighbourhood(int e) {
        std::vector<int>fw;
        double threshold = eps + points_sorted[e].distance;
        for (int i = e + 1; i < points_sorted.size(); i++) {
            // if q.dist > forwardThreshold then 
            // break;
            if (points_sorted[i].distance > threshold) {
                break;
            }
            // if Distance(q, p) <= Eps then
            // append q to neighbours;
            double dist = points_sorted[e].getDistanceN(points_sorted[i], dimensions);
            //increment distance calculation number
            points_sorted[e].noDistances++;
            if (dist <= (eps)) {
                fw.push_back(i);
            }
        }

        return fw;
    }
    vector<int> TI_bwNeighbourhood(int e) {
        ///backward neighbourhood
        vector<int>bw;
        double threshold = eps - points_sorted[e].distance;
        for (int i = e - 1; i >= 0; i--) {
            // if q.dist > forwardThreshold then 
            // break;
            if (points_sorted[i].distance < threshold) {
                break;
            }
            // if Distance2(q, p) <= Eps then 
            // append q to seeds;
            points_sorted[e].noDistances++;
            if (points_sorted[e].getDistanceN(points_sorted[i], dimensions) <= (eps)) {
                bw.push_back(i);
            }
        }
        return bw;
    }
    // check if neighbourhood satisfies the minpoins requirement
    bool isCore(int index) {
        //points_sorted[index].type = 1;
        return points_sorted[index].noOfNeighbours >= minPoints;
    }
    bool isBorder(int ind) {
        for (size_t i = 1; i < neighbours[(int)ind].size(); i++) {
            if (isCore(neighbours[(int)ind][(int)i])) {
                return true;
            }
        }
    }
};
//class reading datapoints from files
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
        string coordinates;

        while (!fin.eof()) {
            fin >> label >> coordinates;
            vector<double> elements;
            std::stringstream ss(coordinates);
            while (ss.good())
            {

                string substr;
                getline(ss, substr, ',');

                elements.push_back(atof(substr.c_str()));
            }
            points.push_back({ index,elements,1, NOT_CLASSIFIED,label });
            index++;
        }

    }
    vector<Point> getPoints() {
        return points;
    }
};


int main(int argc, char* argv[])
{
    auto tstart = chrono::high_resolution_clock::now();

    //path, epsilon and minp are command line arguments
    std::string path = argv[1];
    //std::string path = "simple";
    InputReader input = InputReader(path);
    auto wend = chrono::high_resolution_clock::now();
    auto read_time = chrono::duration_cast<chrono::microseconds>(wend - tstart).count();
    cout << "Read time: " << read_time << " microseconds" << endl;

    double epsilon = atof(argv[2]);
    int minp = atoi(argv[3]);
    TiDbScan alg = TiDbScan(epsilon, minp, input.getPoints());

    auto rstart = chrono::high_resolution_clock::now();
    alg.run();
    auto stop = chrono::high_resolution_clock::now();
    auto time_run = chrono::duration_cast<chrono::microseconds>(stop - rstart).count();
    cout << "Run time (main): " << time_run << " microseconds" << endl;
    auto time_whole = chrono::duration_cast<chrono::microseconds>(stop - tstart).count();
    cout << "Exec time: " << time_whole;
}

// Uruchomienie programu: Ctrl + F5 lub menu Debugowanie > Uruchom bez debugowania
// Debugowanie programu: F5 lub menu Debugowanie > Rozpocznij debugowanie

