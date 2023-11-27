#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include <unordered_set>
#include <limits>

#define COLS 2

using namespace std;

struct Point {
    double x, y;
};

//calculate euclidean distance
double calculateDistance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to create a graph based on the input coordinates
vector<vector<double>> createMatrix(const vector<Point>& points, int size) {
    int n = size;
    vector<vector<double>> matrix(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Calculate Euclidean distance between points i and j
                matrix[i][j] = calculateDistance(points[i], points[j]);
            }
        }
    }

    // Print the distance matrix
   /*  cout << "Distance Matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    } */

    return matrix;
}

vector<int> greedyTSP(vector<vector<double>> matrix, int numCities) {

    vector<int> tour;
    vector<bool> visited(numCities, false);

    int start = 0;
    tour.push_back(start);
    visited[start] = true;

    while(tour.size() < numCities) {
        int nearestCity = -1;
        double minDistance = numeric_limits<double>::max();

        for(int i = 0; i < numCities; i++) {
            if(!visited[i] && matrix[start][i] < minDistance) {
                minDistance = matrix[start][i];
                nearestCity = i;
            }
        }

        if(nearestCity != -1) {
            tour.push_back(nearestCity);
            visited[nearestCity] = true;
            start = nearestCity;
        }
    }

    return tour;
}

int main() {

    int numberOfPoints;

    cin >> numberOfPoints;
    vector<Point> coordinates;

    for (int i = 0; i < numberOfPoints; i++) {
        Point point;
        std::cin >> point.x;
        std::cin >> point.y;
        //add every coordinate to the end
        coordinates.push_back(point);
    }

    vector<vector<double>> matrix = createMatrix(coordinates, numberOfPoints);
    vector<int> tour = greedyTSP(matrix, numberOfPoints);

    for (int j = 0; j < tour.size(); ++j) {
         cout << tour[j] << endl;
    }
    //cout << endl;
 

    return 0;
}
