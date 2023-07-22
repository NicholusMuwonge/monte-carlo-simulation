#include <stdio.h>
#include <omp.h>
#include "monte_carlo.h"

int main() {
    int global_number_of_simulation_points = 500;
    int num_threads_list[] = {1, 2, 4, 8};

    printf("Function: y = x ** x\n");

    printf("\nSequential Code:\n");
    double start_time = omp_get_wtime();
    double result_seq = monte_carlo_integration(global_number_of_simulation_points, 1); 
    double end_time = omp_get_wtime();
    printf("Approximate integral value (Sequential): %lf\n", result_seq);
    printf("Execution Time (Sequential): %lf seconds\n", end_time - start_time);

    printf("\nParallel Code:\n");
    for (int i = 0; i < sizeof(num_threads_list) / sizeof(num_threads_list[0]); i++) {
        int num_threads = num_threads_list[i];
        printf("Number of threads loaded: %d\n", num_threads);

        start_time = omp_get_wtime();
        double result_par = monte_carlo_integration(global_number_of_simulation_points, num_threads);
        end_time = omp_get_wtime();

        printf("Approximate integral value (Parallel, %d Threads): %lf\n", num_threads, result_par);
        printf("Execution Time (Parallel, %d Threads): %lf seconds\n", num_threads, end_time - start_time);

        double speedup = (end_time - start_time) / (end_time - start_time);
        double efficiency = speedup / num_threads;
        printf("Speedup: %lf\n", speedup);
        printf("Efficiency: %lf\n", efficiency);
        printf("\n");
    }

    printf("\nfunction: y = 2x\n");

    printf("\nSequential Code:\n");
    start_time = omp_get_wtime();
    result_seq = monte_carlo_integration(global_number_of_simulation_points, 1);
    end_time = omp_get_wtime();
    printf("Approximate integral value (Sequential): %lf\n", result_seq);
    printf("Execution Time (Sequential): %lf seconds\n", end_time - start_time);

    // Parallel code performance
    printf("\nParallel Code:\n");
    for (int i = 0; i < sizeof(num_threads_list) / sizeof(num_threads_list[0]); i++) {
        int num_threads = num_threads_list[i];
        printf("Number of threads the task is distributed amongst: %d\n", num_threads);

        start_time = omp_get_wtime();
        double result_par = monte_carlo_integration(global_number_of_simulation_points, num_threads);
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


    printf("\nApproximating Pi(π) value using Monte Carlo simulation:\n");
    int no_of_simulation_points = 10000;
    int num_threads_pi = 8;

    double start_time_pi = omp_get_wtime();
    double cumulative_sum_from_all_threads = monte_carlo_pi(no_of_simulation_points, num_threads_pi);
    double end_time_pi = omp_get_wtime();

    printf("Simulation generated Pi(π) value of: %lf\n", cumulative_sum_from_all_threads);
    printf("Time taken to generate Pi(π) value: %lf seconds\n", end_time_pi - start_time_pi);

    return 0;
}
