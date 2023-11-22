// Code is based on pseudocode from this website: https://www2.seas.gwu.edu/~simhaweb/champalg/tsp/tsp.html

// test1 good output based om best savings
//7
// 1
// 9
// 6
// 3
// 2
// 8
// 5
// 0
// 4


// g++ -std=c++11 -o tsp tsp.cpp
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include <unordered_set>


using namespace std;

struct Point {
    double x, y;
};

// Function to compute savings between vertices i and j
double computeSavings(const vector<vector<double>>& graph, int i, int j, int hub) {
    // cout << "i: " << i << endl;
    // cout << "j: " << j << endl;
    // cout << "Savings: " << graph[i][hub] + graph[j][hub] - graph[i][j]<< endl;
    return graph[i][hub] + graph[j][hub] - graph[i][j];

}

// Function to check if adding the segment (i, j) creates a cycle
bool createsCycle(const unordered_set<int>& vertices, int i, int j) {
    return vertices.count(i) != 0 && vertices.count(j) != 0;
}
/// Function to perform the savings-based heuristic
void savingsHeuristic(const vector<vector<double>>& graph, int hub) {
    int n = graph.size();

    // Step 2: VH = V - {h}
    unordered_set<int> VH;
    for (int v = 1; v < n; ++v) {
        VH.insert(v);
    }

    unordered_map<int, int> degreeCount;

    // Step 3-6: Compute savings and sort in increasing order
    vector<pair<double, pair<int, int>>> savingsList;
    for (int i = 1; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {  // Start from i + 1 to avoid duplicates
            double savingsValue = computeSavings(graph, i, j, hub);
            savingsList.push_back({savingsValue, {i, j}});
    }
}

    sort(savingsList.begin(), savingsList.end());



    // Step 7-19: Perform savings-based heuristic
    vector<pair<int, int>> partialTour;
    unordered_set<int> connectedComponents;

    while (!savingsList.empty()) {
        const auto& savingsPair = savingsList.back();
        savingsList.pop_back();

        int i = savingsPair.second.first;
        int j = savingsPair.second.second;

            // cout << "i: " << i << endl;
            // cout << "j: " << j << endl;

            // print connectedComponents
            // cout << "connected comps" << endl;  
            // for (auto& it : connectedComponents) {
            //     cout << it << endl;
            // }

        // Step 9: Check if shortcut does not create a cycle and degree(v) <= 2 for all v
        if (!createsCycle(connectedComponents, i, j)) {
            // Step 10-16: Add segment to partial tour and update VH and connectedComponents
                        // cout << "VH.size(): " << VH.size() << endl;

            // cout << "included"<< endl;
            // cout << "i: " << i << endl;
            // cout << "j: " << j << endl;

            // print degreeCount
            // for (auto& it : degreeCount) {
            //     cout << it.first << " " << it.second << endl;
            // }
            
             // Step 10-16: Add segment to partial tour and update VH and connectedComponents

            // Step 10: Check if degree(v) >= 2
            if (degreeCount[i] < 2 && degreeCount[j] < 2) {
                partialTour.push_back({i, j});
                // connectedComponents.insert(i);
                // connectedComponents.insert(j);

                // Update degreeCount
                degreeCount[i]++;
                degreeCount[j]++;
            }

            if (degreeCount[i] == 2) {
                VH.erase(i);
            }

            if (degreeCount[j] == 2) {
                VH.erase(j);
            }

        }
    }
    cout << "test" << endl; 
    // print degreeCount
    for (auto& it : degreeCount) {
        if (it.second < 2) {
            partialTour.push_back({hub, it.first});
        }
        // cout << it.first << " " << it.second << endl;
    }

    // Step 19: Stitch together remaining two vertices and hub into final tour
    // for (int v : VH) {
    //     partialTour.push_back({hub, v});
    // }

    // Output the results
    double total = 0;
    for (const auto& segment : partialTour) {
        total += graph[segment.first][segment.second];
        // cout << segment.first << endl;
        // cout << segment.second << endl;

        cout << segment.first << " " << segment.second << endl;
    }
    cout << "Total: " << total << endl;
}

// Function to calculate Euclidean distance between two points
double calculateDistance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to create a graph based on the input coordinates
vector<vector<double>> createGraph(const vector<Point>& points, int size) {
    int n = size;
    vector<vector<double>> graph(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Calculate Euclidean distance between points i and j
                graph[i][j] = calculateDistance(points[i], points[j]);
            }
        }
    }

    // Print the distance matrix
// cout << "Distance Matrix:" << endl;
// for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < n; ++j) {
//         cout << graph[i][j] << " ";
//     }
//     cout << endl;
// }

    return graph;
}

int main(int argc, char** argv) {
// int main() {

    int size;

    std::cin >> size;
    std::vector<Point> coordinates;

    for (int i = 0; i < size; i++) {
        Point point;
        std::cin >> point.x;
        std::cin >> point.y;
        coordinates.push_back(point);
    }

    // int size = 10;
    // int size2 = 6;
    // vector<Point> coordinates = {
    //     {95.0129, 61.5432},
    //     {23.1139, 79.1937},
    //     {60.6843, 92.1813},
    //     {48.5982, 73.8207},
    //     {89.1299, 17.6266},
    //     {76.2097, 40.5706},
    //     {45.6468, 93.5470},
    //     {1.8504, 91.6904},
    //     {82.1407, 41.0270},
    //     {44.4703, 89.3650}
    // };

    // vector<Point> coordinates2 = {
    //     {0.0000, 2.0000},
    //     {4.0000, 2.0000},
    //     {8.0000, 2.0000},
    //     {0.0000, 4.0000},
    //     {4.0000, 0.0000},
    //     {8.0000, 4.0000}
    // };


    // Corner case, test 6 ==> Only one point
    if (size == 1){
        std::cout << 0 << std::endl;
        return 0;
    }

    if (size == 2){
        std::cout << 0 << std::endl;
        std::cout << 1 << std::endl;
        return 0;
    }

    if (size == 3){
        std::cout << 0 << std::endl;
        std::cout << 1 << std::endl;
        std::cout << 2 << std::endl;
        return 0;
    }


    vector<vector<double>> graph = createGraph(coordinates, size);

    // vector<vector<double>> graph = createGraph(coordinates2, size2);


    int hub = 0;

    // Perform savings-based heuristic
    savingsHeuristic(graph, hub);

    return 0;

}