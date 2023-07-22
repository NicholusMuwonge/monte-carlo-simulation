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
int within_circle_radius(double x, double y) {
    return (x * x + y * y) <= 1.0;
}

// Function to check if a point (x, y) lies inside the shaded region under the curve y = f(x)
int is_inside_region(double x, double y) {
    return (y <= f(x));
}

// Multiplicative Linear Congruential Generator (MLCG)
double mlcg(double seed, double lambda, double P) {
    return fmod(lambda * seed, P);
}

// Monte Carlo simulation for numerical integration of the function y = f(x) over the interval [0, 10]
double monte_carlo_integration(int number_of_simulations, int num_threads) {
    int inside_shaded_region_count = 0;
    double domain_min_x = 0.0;
    double domain_max_x = 10.0;
    double range_min_y = 0.0;
    double range_max_y = 10.0;

    int seed = time(NULL); // Initialize the seed for MLCG

    // Set the number of threads for the parallel region
    omp_set_num_threads(num_threads);

    #pragma omp parallel for reduction(+:inside_shaded_region_count)
    for (int i = 0; i < number_of_simulations; i++) {
        // Generate random (x, y) coordinates within the specified domain and range using MLCG
        seed = mlcg(seed, 3, 2147483647); // Using λ = 3 and P(Modulus) = 2^31 - 1
        double x = domain_min_x + ((double)seed / 2147483647.0) * (domain_max_x - domain_min_x);
        
        seed = mlcg(seed, 3, 2147483647);
        double y = range_min_y + ((double)seed / 2147483647.0) * (range_max_y - range_min_y);

        // Check if the point (x, y) lies inside the shaded region under the curve y = f(x)
        if (is_inside_region(x, y)) {
            inside_shaded_region_count++; // Increment the count of points inside the region
        }
    }

    // Calculate the ratio of points inside the region to total points generated
    double ratio_inside = (double)inside_shaded_region_count / number_of_simulations;

    // Calculate the total area of the region (domain x range)
    double total_area = (domain_max_x - domain_min_x) * (range_max_y - range_min_y);

    // Approximate the integral of the function y = f(x) using the Monte Carlo simulation
    double approximate_integral = ratio_inside * total_area;
    return approximate_integral;
}

// Monte Carlo simulation to approximate the value of π by generating random points inside a unit square
double monte_carlo_pi(int number_of_simulations, int num_threads) {
    int count_inside_circle = 0;
    double domain_min_x = -1.0;
    double domain_max_x = 1.0;
    double range_min_y = -1.0;
    double range_max_y = 1.0;

    double seed = time(NULL); // Initialize the seed for MLCG

    omp_set_num_threads(num_threads);

    #pragma omp parallel for reduction(+:count_inside_circle)
    for (int i = 0; i < number_of_simulations; i++) {
        // Generate random (x, y) coordinates within the specified domain and range using MLCG
        seed = mlcg(seed, 3, 2147483647); // Using λ = 3 and P = 2^31 - 1
        double x = domain_min_x + (seed / 2147483647.0) * (domain_max_x - domain_min_x);
        
        seed = mlcg(seed, 3, 2147483647);
        double y = range_min_y + (seed / 2147483647.0) * (range_max_y - range_min_y);

        // Check if the point (x, y) lies inside the unit circle
        if (within_circle_radius(x, y)) {
            count_inside_circle++; 
        }
    }

    // Calculate the ratio of points inside the unit circle in relation to the total points generated
    double ratio_inside_circle = (double)count_inside_circle / number_of_simulations;

    // Approximate the value of Pi(π) using the Monte Carlo simulation
    double generated_approximated_pi = 4.0 * ratio_inside_circle;
    return generated_approximated_pi;
}
