#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define ARRAY_SIZE 10 // Adjust the size of the array as needed
#define RANDOM_NUMBER_MAX_SIZE ARRAY_SIZE * 100 // Max random number size that is put in the arrays

void generateArray(int *array, int start, int end, int arraySize){

    // Generate random numbers within the range and store them in the array
    for (int i = 0; i < arraySize; i++) {
        array[i] = rand() % (end - start + 1) + start;
    }

    // Sort the array in ascending order
    for (int i = 0; i < arraySize - 1; i++) {
        for (int j = 0; j < arraySize - i - 1; j++) {
            if (array[j] > array[j + 1]) {
                int temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}

void calculateRange(int total_elements, int num_processes, int rank, int *start_range, int *end_range) {
    
    // Calculate the size of each chunk and the remainder    
    int chunk_size = total_elements / num_processes;
    int remainder = total_elements % num_processes;

    // Initialize the starting point of the range
    int start = total_elements - 1;

    // Adjust the starting point based on the rank of the process
    for (int i = 0; i < rank; i++) {

        // Move the start backward by chunk size
        start -= chunk_size;
        if (remainder > 0) {

            // Adjust for any remaining elements
            start--;
            remainder--;
        }
    }

    // Calculate the end point of the range
    int end = start - chunk_size + 1;
    if (remainder > 0) {
        end--;
    }

    // Assign the calculated start and end points to the output variables
    *start_range = end;
    *end_range = start;

}

void printArray(int *array, int size, int rank) {

    // Print the rank of the process and indicate it's an array
    printf("Rank: %d Array: ", rank);

    // Begin printing the array representation
    printf("[ ");

    // Iterate through the elements of the array
    for (int i = 0; i < size; i++) {

        // Print the current element
        printf("%d ", array[i]);

        // Insert a newline character every 10 elements for better readability
        if ((i + 1) % 10 == 0)
            printf("\n  ");
    }

    // End the array representation
    printf("]\n");
}

void addArrayToEndOfArray(int *array1, int arraySize1, int *array2, int arraySize2, int rank, int num_processes){
    int start, end;

    // Calculate the start and end indices
    start = (num_processes - rank - 1) * ARRAY_SIZE;
    end = start + ARRAY_SIZE - 1;

    // Add the new elements to the new array
    for(int i = start ; i <= end ; i++ ){
        
        // Add new values to array
        array1[i] = array2[i - start];
    }
}

// void addArrayToEndOfArray(int *array1, int arraySize1, int *array2, int arraySize2, int rank, int num_processes) {
//     int start, end;
//     int elements_per_process = arraySize1 / num_processes;
//     int remainder = arraySize1 % num_processes;

//     start = rank * elements_per_process + (rank < remainder ? rank : remainder);
//     end = start + elements_per_process - 1 + (rank < remainder ? 1 : 0);

//     // Copy elements from array2 to array1
//     for (int i = 0; i < arraySize2; i++) {
//         array1[start + i] = array2[i];
//     }
// }

int main(int argc, char *argv[]) {

    // Seed the random number generator once to ensure randomness
    srand(time(NULL));

    int rank, size;
    double start_time, end_time, elapsed_time, tick;
    int *g_sortedIntegerArray1;
    int *g_sortedIntegerArray2;

    // Define variables to hold the size of the global arrays and set their value
    int size_g_sortedIntegerArray1;
    int size_g_sortedIntegerArray2;
        
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Returns the precision of the results returned by MPI_Wtime
    tick = MPI_Wtick();

    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    int start_range, end_range;

    // Calculate the range of elements for the current process
    calculateRange(RANDOM_NUMBER_MAX_SIZE, size, rank, &start_range, &end_range);

    // Initialize 2 arrays to hold the random integers
    int *sortedRandomIntegerArray1 = (int *)malloc(ARRAY_SIZE * sizeof(int));
    int *sortedRandomIntegerArray2 = (int *)malloc(ARRAY_SIZE * sizeof(int));

    // Generate 2 arrays of random integers for the current process
    generateArray(sortedRandomIntegerArray1, start_range, end_range, ARRAY_SIZE);
    generateArray(sortedRandomIntegerArray2, start_range, end_range, ARRAY_SIZE);

    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    end_time = MPI_Wtime(); 

    MPI_Status status;

    // Calculate the size of the global arrays
    size_g_sortedIntegerArray1 = ARRAY_SIZE * size;
    size_g_sortedIntegerArray2 = ARRAY_SIZE * size;

    // Allocate memory for global sorted array
    g_sortedIntegerArray1 = (int *)malloc(size_g_sortedIntegerArray1 * sizeof(int));
    g_sortedIntegerArray2 = (int *)malloc(size_g_sortedIntegerArray2 * sizeof(int));

    // if rank is 0, collect all data from all other processes
    if (rank == 0) {

        // Array to store received arrays
        int receivedArrays1[size][ARRAY_SIZE]; 
        int receivedArrays2[size][ARRAY_SIZE]; 
        
        // Receive arrays from all processes except for rank 0
        for (int source = size - 1; source > 0; source--) {

            MPI_Recv(receivedArrays1[source], ARRAY_SIZE, MPI_INT, source, 1, MPI_COMM_WORLD, &status);

            // Add received array to the end the of the array
            addArrayToEndOfArray(g_sortedIntegerArray1, size_g_sortedIntegerArray1, receivedArrays1[source], ARRAY_SIZE, source, size);

            MPI_Recv(receivedArrays2[source], ARRAY_SIZE, MPI_INT, source, 2, MPI_COMM_WORLD, &status);

            // Add received array to the end the of the array
            addArrayToEndOfArray(g_sortedIntegerArray2, size_g_sortedIntegerArray2, receivedArrays2[source], ARRAY_SIZE, source, size);

        }

        // Add array from rank 0 to the end the of the array
        addArrayToEndOfArray(g_sortedIntegerArray1, (ARRAY_SIZE * size), sortedRandomIntegerArray1, ARRAY_SIZE, rank, size);

        // Add array from rank 0 to the end the of the array
        addArrayToEndOfArray(g_sortedIntegerArray2, (ARRAY_SIZE * size), sortedRandomIntegerArray2, ARRAY_SIZE, rank, size);

        printf("********** Global Sorted Array 1 **********\n");
        printArray(g_sortedIntegerArray1, ARRAY_SIZE * size, rank);

        printf("********** Global Sorted Array 2 **********\n");
        printArray(g_sortedIntegerArray2, ARRAY_SIZE * size, rank);

        for (int i = 0; i < size; i ++){

            MPI_Send(g_sortedIntegerArray1, size_g_sortedIntegerArray1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(g_sortedIntegerArray2, size_g_sortedIntegerArray2, MPI_INT, i, 2, MPI_COMM_WORLD);
        }

    
    } else {
    
        // Send the two arrays generated by other processes to rank 0 with tag different tags
        MPI_Send(sortedRandomIntegerArray1, ARRAY_SIZE, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(sortedRandomIntegerArray2, ARRAY_SIZE, MPI_INT, 0, 2, MPI_COMM_WORLD);

        MPI_Recv(g_sortedIntegerArray1, size_g_sortedIntegerArray1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        MPI_Recv(g_sortedIntegerArray2, size_g_sortedIntegerArray2, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);


    }
    
    // Calculate the elapsed time
    elapsed_time = end_time - start_time;

    // if rank is 0, print out the time analysis
    if (rank == 0) {
        printf("\n********** Array Creation Time **********\n");
        // printf("Size of arrays to be merged: %d\n", (ARRAY_SIZE * size));
        printf("Total processes: %d\n", size);
        printf("Total computation time: %e seconds\n", elapsed_time);
        printf("Computation time per process: %e seconds\n", elapsed_time / size);
        printf("Resolution of MPI_Wtime: %e seconds\n\n", tick);
    }

    // Flush all print statements before next section of code
    fflush(stdout);

    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    // MERGE ARRAYS

    // split up first array evenly

    // int remainder = arraySize1 % size;
    // int quotient = arraySize1 / size;

    // Calculate the range of elements for the current process
    calculateRange(size_g_sortedIntegerArray1, size, rank, &start_range, &end_range);

    // printf("Rank: %d size_g_sortedIntegerArray1: %d size_g_sortedIntegerArray2: %d\n", rank, size_g_sortedIntegerArray1, size_g_sortedIntegerArray2);

    printf("Rank: %d Start: %d End: %d\n", rank, start_range, end_range);

    // Flush all print statements before next section of code
    fflush(stdout);
    
    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    end_time = MPI_Wtime();

    // Calculate the elapsed time
    elapsed_time = end_time - start_time;

    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    // if rank is 0, print out the time analysis
    if (rank == 0) {
        printf("\n********** Array Merge Time **********\n");
        printf("Total processes: %d\n", size);
        printf("Total computation time: %e seconds\n", elapsed_time);
        printf("Computation time per process: %e seconds\n", elapsed_time / size);
        printf("Resolution of MPI_Wtime: %e seconds\n", tick);
    }

    // Clean up
    free(sortedRandomIntegerArray1);
    free(sortedRandomIntegerArray2);
    free(g_sortedIntegerArray1);
    free(g_sortedIntegerArray2);

    MPI_Finalize();

    return 0;
}
