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

void calculateRangeBasedOnRankAndSize(int total_elements, int num_processes, int rank, int *start_range, int *end_range) {
    
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

void calculateRangeBasedOnArray(int *array1, int *array2, int arraySize2, int array1_start_range, int array1_end_range, int *array2_start_range, int *array2_end_range){
    
    int start = -1; // Initialize start index
    int end = -1;   // Initialize end index

    int upper_value = array1[array1_end_range];
    

    // Assign the calculated start and end points to the output variables
    if (array1_start_range == 0){

        *array2_start_range = 0;

    }else {

        // Set the lower value to one index less to account for values that could be between
        int lower_value = array1[(array1_start_range - 1)];  

        // Find the leftmost index with value greater than lower_value
        for (int i = 0; i < arraySize2; i++) {
            if (array2[i] > lower_value) {
                start = i;
                break;
            }
        }

        *array2_start_range = start;
    }

    // Find the rightmost index with value less than or equal to upper_value
    for (int i = arraySize2 - 1; i >= 0; i--) {
        if (array2[i] <= upper_value) {
            end = i;
            break;
        }
    }

    
    *array2_end_range = end;
}

void mergeTwoSortedArrays(int *sortedArray1, int *sortedArray2, int array1_start_range, int array1_end_range, int array2_start_range, int array2_end_range, int *mergedArray){
        
    int i = array1_start_range, j = array2_start_range, k = 0;

    // Merge elements from sortedArray1 and sortedArray2 into mergedArray
    while (i <= array1_end_range && j <= array2_end_range) {
        if (sortedArray1[i] <= sortedArray2[j]) {
            mergedArray[k] = sortedArray1[i];
            k++;
            i++;
        } else {
            mergedArray[k] = sortedArray2[j];
            k++;
            j++;
        }
    }

    // Copy the remaining elements of sortedArray1, if any
    while (i <= array1_end_range) {
        mergedArray[k] = sortedArray1[i];
        k++;
        i++;
    }

    // Copy the remaining elements of sortedArray2, if any
    while (j <= array2_end_range) {
        mergedArray[k] = sortedArray2[j];
        k++;
        j++;
    }

}

int main(int argc, char *argv[]) {

    // Seed the random number generator once to ensure randomness
    srand(time(NULL));

    // initialize random sorted array and their lengths
    int n_online_random_array1[] = {69, 179, 227, 249, 254, 398, 556, 570, 619, 776, 789, 813, 862, 927, 988, 998, 1044, 1197, 1657, 1699, 1820, 1973, 1974, 2044, 2051, 2163, 2392, 2398, 2580, 2621, 2751, 2755, 2955, 3125, 3199, 3272, 3303, 3490, 3540, 3551, 3606, 3616, 4021, 4066, 4206, 4214, 4252, 4380, 4624, 4751, 4788, 5004, 5078, 5114, 5459, 5593, 5839, 5930, 5992, 5994, 6005, 6011, 6063, 6127, 6151, 6165, 6298, 6317, 6411, 6553, 6599, 6630, 6779, 7455, 7521, 7571, 7652, 7704, 7776, 7868, 7869, 8322, 8454, 8475, 8507, 8530, 8538, 8546, 8608, 8711, 8824, 8910, 8941, 9107, 9195, 9206, 9404, 9504, 9505, 9723};
    int n_online_random_array2[] = {32, 59, 120, 166, 358, 491, 539, 774, 787, 914, 943, 1108, 1222, 1234, 1249, 1305, 1309, 1344, 1526, 1575, 1772, 1940, 1964, 2185, 2197, 2322, 2411, 2412, 2486, 2625, 2965, 2972, 3124, 3126, 3192, 3355, 3391, 3534, 3616, 3632, 3683, 3788, 3849, 3989, 4026, 4111, 4174, 4413, 4513, 4752, 4756, 4795, 4919, 4946, 4989, 5125, 5493, 5541, 5659, 5725, 5754, 5778, 5833, 5837, 5864, 5876, 5987, 6024, 6290, 6313, 6340, 6534, 6572, 6587, 6702, 6759, 6808, 6822, 6840, 7289, 7451, 7527, 7697, 7735, 7831, 7931, 7933, 7967, 8004, 8045, 8113, 8237, 8494, 8516, 8606, 8664, 9020, 9628, 9635, 9851};
    int *online_random_array1 = n_online_random_array1;
    int *online_random_array2 = n_online_random_array2;
    int size_online_random_array1 = 100;
    int size_online_random_array2 = 100;

    // initialize random sorted array that will be generated locally
    int *g_sortedIntegerArray1;
    int *g_sortedIntegerArray2;

    // Define variables to hold the size of the locally generate arrays
    int size_g_sortedIntegerArray1;
    int size_g_sortedIntegerArray2;

    int rank, size;
    double start_time, end_time, elapsed_time, tick;
        
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Returns the precision of the results returned by MPI_Wtime
    tick = MPI_Wtick();

    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    int start_range, end_range;

    // LOCALLY GENERATE RANDOM SORTED INTEGER ARRAYS

    // Calculate the range of elements for the current process
    calculateRangeBasedOnRankAndSize(RANDOM_NUMBER_MAX_SIZE, size, rank, &start_range, &end_range);

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
        addArrayToEndOfArray(g_sortedIntegerArray1, size_g_sortedIntegerArray1, sortedRandomIntegerArray1, ARRAY_SIZE, rank, size);

        // Add array from rank 0 to the end the of the array
        addArrayToEndOfArray(g_sortedIntegerArray2, size_g_sortedIntegerArray2, sortedRandomIntegerArray2, ARRAY_SIZE, rank, size);

        // Testing purposes only
        printf("********** Global Sorted Array 1 **********\n");
        printArray(g_sortedIntegerArray1, size_g_sortedIntegerArray1, rank);

        // Testing purposes only
        printf("********** Global Sorted Array 2 **********\n");
        printArray(g_sortedIntegerArray2, size_g_sortedIntegerArray2, rank);

        // Send the global sorted arrays to all other processes
        for (int i = 0; i < size; i ++){

            MPI_Send(g_sortedIntegerArray1, size_g_sortedIntegerArray1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(g_sortedIntegerArray2, size_g_sortedIntegerArray2, MPI_INT, i, 2, MPI_COMM_WORLD);
        }

    
    } else {
    
        // Send the two arrays generated by other processes to rank 0 with tag different tags
        MPI_Send(sortedRandomIntegerArray1, ARRAY_SIZE, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(sortedRandomIntegerArray2, ARRAY_SIZE, MPI_INT, 0, 2, MPI_COMM_WORLD);

        // Receive global sorted arrays from process 0
        MPI_Recv(g_sortedIntegerArray1, size_g_sortedIntegerArray1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(g_sortedIntegerArray2, size_g_sortedIntegerArray2, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);

    }
    
    // Calculate the elapsed time
    elapsed_time = end_time - start_time;

    // if rank is 0, print out the time analysis of array creation
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

    // initialize variables to hold start and end of each array that each process will work on
    int array1_start_range, array1_end_range, array2_start_range, array2_end_range;
    
    // ########## Testing with locally generated array ##########

    // // Calculate the range of elements for array 1 in the current process
    // calculateRangeBasedOnRankAndSize(size_g_sortedIntegerArray1, size, rank, &array1_start_range, &array1_end_range);

    // // Calculate the range of elements for array 2 in the current process
    // calculateRangeBasedOnArray(g_sortedIntegerArray1, g_sortedIntegerArray2, size_g_sortedIntegerArray2, array1_start_range, array1_end_range, &array2_start_range, &array2_end_range);

    // Testing purposes only
    // printf("\nLocal:\nRank: %d Array 1 Start: %d Array 1 Start Value: %d Array 1 End: %d Array 1 End Value: %d\nRank: %d Array 2 Start: %d Array 2 Start Value: %d Array 2 End: %d Array 2 End Value: %d\n", rank, array1_start_range, g_sortedIntegerArray1[array1_start_range], array1_end_range, g_sortedIntegerArray1[array1_end_range], rank, array2_start_range, g_sortedIntegerArray2[array2_start_range], array2_end_range, g_sortedIntegerArray2[array2_end_range]);

    // // Flush all print statements before next section of code
    // fflush(stdout);

    // // Merge the two arrays together
    // mergeTwoSortedArrays(g_sortedIntegerArray1, size_g_sortedIntegerArray1, g_sortedIntegerArray2, size_g_sortedIntegerArray2, merged_array);

    // ########## Testing with locally generated array ##########



    // ########## Testing with random online array ##########

    // Calculate the range of elements for array 1 in the current process
    calculateRangeBasedOnRankAndSize(size_online_random_array1, size, rank, &array1_start_range, &array1_end_range);

    // Calculate the range of elements for array 2 in the current process
    calculateRangeBasedOnArray(online_random_array1, online_random_array2, size_online_random_array2, array1_start_range, array1_end_range, &array2_start_range, &array2_end_range);

    // Testing purposes only
    // printf("\nOnline:\nRank: %d Array 1 Start: %d Array 1 Start Value: %d Array 1 End: %d Array 1 End Value: %d\nRank: %d Array 2 Start: %d Array 2 Start Value: %d Array 2 End: %d Array 2 End Value: %d\n", rank, array1_start_range, online_random_array1[array1_start_range], array1_end_range, online_random_array1[array1_end_range], rank, array2_start_range, online_random_array2[array2_start_range], array2_end_range, online_random_array2[array2_end_range]);

    // Flush all print statements before next section of code
    fflush(stdout);

    int *merged_array;
    int size_merged_array = (array1_end_range - array1_start_range) + (array2_end_range - array2_start_range) + 2;

    // Allocate memory for merged array
    merged_array = (int *)malloc(size_merged_array * sizeof(int));

    // Merge the two arrays together
    mergeTwoSortedArrays(online_random_array1, online_random_array2, array1_start_range, array1_end_range, array2_start_range, array2_end_range, merged_array);
    
    // ########## Testing with random online array ##########


    // if rank is 0, collect all merged arrays from all other processes
    if (rank == 0) {

        // initialize final sorted array
        int *final_sortedIntegerArray;

        // Calculate the size of the final sorted array 
        int size_final_sortedIntegerArray = (size_g_sortedIntegerArray1 + size_g_sortedIntegerArray2);

        // Allocate memory for global sorted array
        final_sortedIntegerArray = (int *)malloc(size_final_sortedIntegerArray * sizeof(int));

        // initialize variables to receive data from other proccess
        int *receivedMergedArray;
        int receivedMergedArraySize;

        // variable to help with inserting into final sorted array
        int finalArrayCurrPostion = 0;

        // Receive arrays from all processes except for rank 0 from largest to smallest
        for (int source = size - 1; source > 0; source--) {
            
            // receive merged array size
            MPI_Recv(&receivedMergedArraySize, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);

            // allocate space for array
            receivedMergedArray = (int *)malloc(receivedMergedArraySize * sizeof(int));

            // receive merged array
            MPI_Recv(receivedMergedArray, receivedMergedArraySize, MPI_INT, source, 1, MPI_COMM_WORLD, &status);

            // Testing purposes only
            printf("receivedMergedArraySize: %d\n", receivedMergedArraySize);
            printArray(receivedMergedArray, receivedMergedArraySize, source);
            

            // Add array to end of array bases on size and finalArrayCurrPostion
        
        }

        // Add array from rank 0 to end of array bases on size and finalArrayCurrPostion
        

        // printf("********** Final Sorted Array 1 **********\n");
        // printArray(final_sortedIntegerArray, size_final_sortedIntegerArray, rank);
    
    } else {

        // Send merged array size to process 0
        MPI_Send(&size_merged_array, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

        // Send merged array to process 0
        MPI_Send(merged_array, size_merged_array, MPI_INT, 0, 1, MPI_COMM_WORLD);

    }

    // Flush all print statements before next section of code
    fflush(stdout);
    
    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    end_time = MPI_Wtime();

    // Calculate the elapsed time
    elapsed_time = end_time - start_time;

    // Ensures all processes will enter the measured section of the code at the same time
    MPI_Barrier(MPI_COMM_WORLD);

    // if rank is 0, print out the time analysis for merging arrays
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
    free(receivedMergedArray);

    MPI_Finalize();

    return 0;
}
