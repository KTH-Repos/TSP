#include <stdio.h>

#define COLS 2

int main() {

    int numberOfPoints;

    scanf("%d", &numberOfPoints);

    float coordinates[numberOfPoints][COLS];

    for (int i = 0; i < numberOfPoints; i++) {
        for (int j = 0; j < COLS; j++) {
            scanf("%f", &coordinates[i][j]);
        }
    }

    /* printf("here is the print\n");
    for (int i = 0; i < numberOfPoints; i++) {
        for (int j = 0; j < COLS; j++) {
            printf("%f\t", coordinates[i][j]);
        }
        printf("\n");
    } */

    return 0;
}
