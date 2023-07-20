#include <stdio.h>
#include <omp.h>
#include "monte_carlo.h"

int main() {
    int N = 500; // Number of random points
    int num_threads_list[] = {1, 2, 4, 8}; // List of number of threads to test

    printf("Function: y = x^2\n");

    // Sequential code performance
    printf("\nSequential Code:\n");
    double start_time = omp_get_wtime();
    double result_seq = monte_carlo_integration(N, 1); // Single thread for sequential code
    double end_time = omp_get_wtime();
    printf("Approximate integral value (Sequential): %lf\n", result_seq);
    printf("Execution Time (Sequential): %lf seconds\n", end_time - start_time);

    // Parallel code performance
    printf("\nParallel Code:\n");
    for (int i = 0; i < sizeof(num_threads_list) / sizeof(num_threads_list[0]); i++) {
        int num_threads = num_threads_list[i];
        printf("Number of Threads: %d\n", num_threads);

        start_time = omp_get_wtime();
        double result_par = monte_carlo_integration(N, num_threads);
        end_time = omp_get_wtime();

        printf("Approximate integral value (Parallel, %d Threads): %lf\n", num_threads, result_par);
        printf("Execution Time (Parallel, %d Threads): %lf seconds\n", num_threads, end_time - start_time);

        // Performance comparison
        double speedup = (end_time - start_time) / (end_time - start_time);
        double efficiency = speedup / num_threads;
        printf("Speedup: %lf\n", speedup);
        printf("Efficiency: %lf\n", efficiency);
        printf("\n");
    }

    // Now, let's use the new function y = 2x for the Monte Carlo simulation
    printf("\nFunction: y = 2x\n");

    // Sequential code performance
    printf("\nSequential Code:\n");
    start_time = omp_get_wtime();
    result_seq = monte_carlo_integration(N, 1); // Single thread for sequential code
    end_time = omp_get_wtime();
    printf("Approximate integral value (Sequential): %lf\n", result_seq);
    printf("Execution Time (Sequential): %lf seconds\n", end_time - start_time);

    // Parallel code performance
    printf("\nParallel Code:\n");
    for (int i = 0; i < sizeof(num_threads_list) / sizeof(num_threads_list[0]); i++) {
        int num_threads = num_threads_list[i];
        printf("Number of Threads: %d\n", num_threads);

        start_time = omp_get_wtime();
        double result_par = monte_carlo_integration(N, num_threads);
        end_time = omp_get_wtime();

        printf("Approximate integral value (Parallel, %d Threads): %lf\n", num_threads, result_par);
        printf("Execution Time (Parallel, %d Threads): %lf seconds\n", num_threads, end_time - start_time);

        // Performance comparison
        double speedup = (end_time - start_time) / (end_time - start_time);
        double efficiency = speedup / num_threads;
        printf("Speedup: %lf\n", speedup);
        printf("Efficiency: %lf\n", efficiency);
        printf("\n");
    }


    printf("\nApproximating π using Monte Carlo simulation:\n");
    int N_pi = 1000000; // Number of random points for π approximation
    int num_threads_pi = 4; // Number of threads for π approximation

    double start_time_pi = omp_get_wtime();
    double approximate_pi = monte_carlo_pi(N_pi, num_threads_pi); // Call the function with correct declaration
    double end_time_pi = omp_get_wtime();

    printf("Approximate value of π: %lf\n", approximate_pi);
    printf("Execution Time (π approximation): %lf seconds\n", end_time_pi - start_time_pi);


    return 0;
}
