#include "monte_carlo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

// Function f(x) = x^2
double f(double x) {
    return x * x;
}

// Function to check if a point (x, y) lies inside the unit circle
int is_inside_circle(double x, double y) {
    return (x * x + y * y) <= 1.0;
}

// Function to check if a point (x, y) lies inside the shaded region under the curve y = f(x)
int is_inside_region(double x, double y) {
    return (y <= f(x));
}

// Monte Carlo simulation for numerical integration of the function y = f(x) over the interval [0, 10]
double monte_carlo_integration(int N, int num_threads) {
    int count_inside = 0;
    double domain_min = 0.0;
    double domain_max = 10.0;
    double range_min = 0.0;
    double range_max = 100.0; // For the new function y = 2x, range is [0, 100]

    srand(time(NULL)); // Seed the random number generator

    // Set the number of threads for the parallel region
    omp_set_num_threads(num_threads);

    // Parallel loop to generate N random points and count the points inside the shaded region
    #pragma omp parallel for reduction(+:count_inside)
    for (int i = 0; i < N; i++) {
        // Generate random (x, y) coordinates within the specified domain and range
        double x = domain_min + (double)rand() / RAND_MAX * (domain_max - domain_min);
        double y = range_min + (double)rand() / RAND_MAX * (range_max - range_min);

        // Check if the point (x, y) lies inside the shaded region under the curve y = f(x)
        if (is_inside_region(x, y)) {
            count_inside++; // Increment the count of points inside the region
        }
    }

    // Calculate the ratio of points inside the region to total points generated
    double ratio_inside = (double)count_inside / N;

    // Calculate the total area of the region (domain x range)
    double total_area = (domain_max - domain_min) * (range_max - range_min);

    // Approximate the integral of the function y = f(x) using the Monte Carlo simulation
    double approximate_integral = ratio_inside * total_area;
    return approximate_integral;
}

// Monte Carlo simulation to approximate the value of π by generating random points inside a unit square
double monte_carlo_pi(int N, int num_threads) {
    int count_inside_circle = 0;
    double domain_min = -1.0;
    double domain_max = 1.0;
    double range_min = -1.0;
    double range_max = 1.0;

    srand(time(NULL)); // Seed the random number generator

    // Set the number of threads for the parallel region
    omp_set_num_threads(num_threads);

    // Parallel loop to generate N random points and count the points inside the unit circle
    #pragma omp parallel for reduction(+:count_inside_circle)
    for (int i = 0; i < N; i++) {
        // Generate random (x, y) coordinates within the specified domain and range
        double x = domain_min + (double)rand() / RAND_MAX * (domain_max - domain_min);
        double y = range_min + (double)rand() / RAND_MAX * (range_max - range_min);

        // Check if the point (x, y) lies inside the unit circle
        if (is_inside_circle(x, y)) {
            count_inside_circle++; // Increment the count of points inside the circle
        }
    }

    // Calculate the ratio of points inside the unit circle to total points generated
    double ratio_inside_circle = (double)count_inside_circle / N;

    // Approximate the value of π using the Monte Carlo simulation
    double approximate_pi = 4.0 * ratio_inside_circle;
    return approximate_pi;
}
