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
#include <stdlib.h>
#include <unordered_set>
#include <unordered_map>
#include <queue>

using namespace std;

struct Point {
    double x, y;
};


// Union find algorithm
class UnionFind {
public:
    unordered_map<int, int> sets;

    // Find the set that our node belongs to and that set's parent
    int find(int node) {

        // If the node is not in the set, add it
        if (sets.find(node) == sets.end()) {
            sets[node] = node;
            return node;
        }

        // If the node is the parent, return parent
        if (sets[node] == node) {
            return node;
        }

        // Recursively iterate through the set to find the parent
        int root = find(sets[node]);

        // Path compression - Create a direct connection from node in set to parent.
        sets[node] = root;
        return root;
    }

    // Add node2 to the set that root belongs to, if they are not already in the same set
    void unionSets(int node1, int node2) {
        int root1 = find(node1);
        int root2 = find(node2);
        
        if (root1 != root2) {
            sets[root1] = root2;
        }
    }
};

// Function to compute savings between vertices i and j
double computeSavings(const vector<vector<double>>& graph, int i, int j, int hub) {
    return graph[i][hub] + graph[j][hub] - graph[i][j];
}


// Check if wanted pair is in our partialTour list
bool checkCycle(vector<pair<int, int>> partialTour, int i, int j, UnionFind& unionFind) {

    int root_i = unionFind.find(i);
    int root_j = unionFind.find(j);
    
    // Return TRUE if they belong in the same set, i.e have the same parent
    return root_i == root_j;
}

void printFinalPath(unordered_map<int, vector<int>> bestPairs) {
    int counter = 0;
    int offset = 0;
    unordered_map<int, int> moveOrder;

    // for(int i=0; i<bestPairs.size(); i = offset){
    //     std::cout << "i: " << i << std::endl; // i should be 5 on the third iteration but we are never going in here after the second iteration.

    //     if(i!=0){
    //         counter++;
    //     }
    //     moveOrder.insert({i, counter});
    //     vector<int> connectedNodes = bestPairs[i];
    //     for (int value: connectedNodes){
    //         if(bestPairs.count(value)!=0){
    //             offset = value;
    //         }
    //     }

    //     bestPairs.erase(i);
    //     std::cout << "offset: " << offset << std::endl;
    //     std::cout << "size of bestPairs: " << bestPairs.size() << std::endl;
    // }

    while (!bestPairs.empty()) {
        // std::cout << "i: " << offset << std::endl;

        if (offset != 0) {
            counter++;
        }

        moveOrder.insert({offset, counter});

        const std::vector<int>& connectedNodes = bestPairs[offset];
        int nextOffset = -1;  // Initialize nextOffset to an invalid value

        for (int value : connectedNodes) {
            if (bestPairs.count(value) != 0) {
                nextOffset = value;
                break;  // Exit the inner loop if the value is found
            }
        }

        // Erase the current offset
        bestPairs.erase(offset);
        offset = nextOffset;
    }

    // TODO: print the moveOrder as Kattis wants it.

    std::vector<std::pair<int, int>> moveOrderVector(moveOrder.begin(), moveOrder.end());
    sort(moveOrderVector.begin(), moveOrderVector.end());

    for (const auto& pair : moveOrderVector) {
        std::cout << pair.second << std::endl;
    }

}

// Function to perform the savings-based heuristic
void savingsHeuristic(const vector<vector<double>>& graph, int hub) {
    int n = graph.size();
    UnionFind unionFind;

    // Kombinera denna med degreecount
    unordered_map<int, vector<int>> connectedComponents;

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

    while (!savingsList.empty()) {
        const auto& savingsPair = savingsList.back();
        savingsList.pop_back();

        int i = savingsPair.second.first;
        int j = savingsPair.second.second;

        //cout << "savings: " << i << " " << j << endl;

        if (!checkCycle(partialTour, i, j, unionFind)) {
        // Step 9: Check if shortcut does not create a cycle and degree(v) <= 2 for all v
        
            
            // Step 10-16: Add segment to partial tour and update VH and connectedComponents

            // Step 10: Check if degree(v) >= 2
            if (degreeCount[i] < 2 && degreeCount[j] < 2) {
                partialTour.push_back({i, j});
                // connectedComponents.insert(i);
                // connectedComponents.insert(j);

                // Add j to the same set as the one i belongs to
                unionFind.unionSets(i, j); 

                connectedComponents[i].push_back(j);
                connectedComponents[j].push_back(i);

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

    for (auto& it : degreeCount) {
        if (it.second < 2) {
            partialTour.push_back({hub, it.first});
            connectedComponents[hub].push_back(it.first);
            connectedComponents[it.first].push_back(hub);
        }
    }

    printFinalPath(connectedComponents);
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

    return graph;
}

int main(int argc, char** argv) {

    int size;

    std::cin >> size;
    std::vector<Point> coordinates;

    for (int i = 0; i < size; i++) {
        Point point;
        std::cin >> point.x;
        std::cin >> point.y;
        coordinates.push_back(point);
    }

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

    int hub = 0;

    // Perform savings-based heuristic
    savingsHeuristic(graph, hub);

    return 0;

}