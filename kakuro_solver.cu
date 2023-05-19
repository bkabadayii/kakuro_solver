#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <stack>
#include <bits/stdc++.h>
#include <math.h>

using namespace std;

enum direction
{
    d_down,
    d_right,
    none
};

#define COORD std::pair<int, int>

// #define DEBUG

int iter = 0;

//////////////////////////////////////////////
// Auxiliary functions for preparing problem //
//////////////////////////////////////////////

void display_arr(int *arr, int n)
{

    cout << "arr: ";

    for (int i = 0; i < n; i++)
    {
        cout << arr[i] << " ";
    }

    cout << endl;
}

void print_coords(COORD start, COORD end)
{

    cout << "Start:" << start.first << "," << start.second << endl;
    cout << "End:" << end.first << "," << end.second << endl;
}

int find_length(COORD start, COORD end, direction dir)
{

    if (dir == d_down)
        return end.first - start.first;
    if (dir == d_right)
        return end.second - start.second;

    return -1;
}

void convert_sol(int **mat, int **&sol_mat, int m, int n)
{

    sol_mat = new int *[m]; // Rows
    for (int i = 0; i < m; i++)
    {
        sol_mat[i] = new int[n]; // Cols
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (mat[i][j] == -2)
                sol_mat[i][j] = -2; // Empty value cell
            else
                sol_mat[i][j] = -1; // Hint or empty cell
        }
    }
}

void print_one_matrix(int **matrix, int m, int n)
{
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++)
    { // rows
        for (int j = 0; j < n; j++)
        { // cols
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

/// Auxiliary functions

struct sum
{
    COORD start;
    COORD end;

    int hint;
    int dir;
    int length;
    int *arr;

    void print_sum()
    {
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << "############################" << endl;
    }

    sum(COORD _start, COORD _end, int _hint, direction _dir) : start(_start), end(_end), hint(_hint), dir(_dir)
    {
        length = find_length(_start, _end, _dir);
        arr = new int[length];
#ifdef DEBUG
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << "############################" << endl;
#endif
    }

    //~sum(){
    // delete arr;
    //}
};

COORD find_end(int **matrix, int m, int n, int i, int j, direction dir)
{ // 0 down 1 right

    if (dir == d_right)
    {
        for (int jj = j + 1; jj < n; jj++)
        {
            if (matrix[i][jj] != -2 || jj == n - 1)
            {
                if (matrix[i][jj] == -2 && jj == n - 1)
                    jj++;
                COORD END = COORD(i, jj);
                return END;
            }
        }
    }

    if (dir == d_down)
    {
        for (int ii = i + 1; ii < m; ii++)
        {
            if (matrix[ii][j] != -2 || ii == m - 1)
            {
                if (matrix[ii][j] == -2 && ii == m - 1)
                    ii++;
                COORD END = COORD(ii, j);
                return END;
            }
        }
    }

    return COORD();
}

vector<sum> get_sums(int **matrix, int m, int n)
{

    vector<sum> sums;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int val = matrix[i][j];
            if (val != -1 && val != -2)
            {
                int hint = val;
                hint = hint / 10;

                if ((hint % 100) == 0)
                {
                    hint = (int)(hint / 100);
                    COORD START = COORD(i, j + 1);
                    COORD END = find_end(matrix, m, n, i, j, d_right);
                    sum _sum = sum(START, END, hint, d_right);
                    sums.push_back(_sum);
                }

                else
                {
                    int div = (int)(hint / 100);
                    int rem = (int)(hint % 100);

                    if (div == 0 && rem != 0)
                    {
                        COORD START = COORD(i + 1, j);
                        COORD END = find_end(matrix, m, n, i, j, d_down);
                        sum _sum = sum(START, END, rem, d_down);
                        sums.push_back(_sum);
                    }

                    if (div != 0 && rem != 0)
                    {
                        COORD START1 = COORD(i + 1, j);
                        COORD START2 = COORD(i, j + 1);
                        COORD END1 = find_end(matrix, m, n, i, j, d_down);
                        COORD END2 = find_end(matrix, m, n, i, j, d_right);
                        sum _sum1 = sum(START1, END1, rem, d_down);
                        sum _sum2 = sum(START2, END2, div, d_right);
                        sums.push_back(_sum1);
                        sums.push_back(_sum2);
                    }
                }
            }
        }
    }
    return sums;
}

void read_matrix(int **&matrix, std::ifstream &afile, int m, int n)
{

    matrix = new int *[m]; // rows

    for (int i = 0; i < m; i++)
    {
        matrix[i] = new int[n]; // cols
    }

    int val;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            afile >> val;
            matrix[i][j] = val;
        }
    }
}

void sol_to_file(int **mat, int **sol_mat, int m, int n)
{

    string fname = "visualize.kakuro";
    ofstream to_write(fname);

    to_write << m << " " << n << "\n";

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat[i][j] != -2)
                to_write << mat[i][j] << " ";
            else
                to_write << sol_mat[i][j] << " ";
        }
        to_write << "\n";
    }

    to_write.close();
}

//////////////////////////////////////////////
// Auxiliary functions for preparing problem //
//////////////////////////////////////////////

///////////////////////////////////////////////////
// Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////

void flatten_sums(vector<sum> sums, int *h_sum_starts_x, int *h_sum_starts_y, int *h_sum_ends_x, int *h_sum_ends_y, int *h_sum_hints, int *h_sum_lengths, int *h_sum_dirs, int no_sums)
{

    for (int i = 0; i < no_sums; i++)
    {

        h_sum_starts_x[i] = sums[i].start.first;
        h_sum_starts_y[i] = sums[i].start.second;

        h_sum_ends_x[i] = sums[i].end.first;
        h_sum_ends_y[i] = sums[i].end.second;

        h_sum_hints[i] = sums[i].hint;
        h_sum_lengths[i] = sums[i].length;

        h_sum_dirs[i] = sums[i].dir;
    }
}

void print_flattened(int *h_sum_starts_x, int *h_sum_starts_y, int *h_sum_ends_x, int *h_sum_ends_y, int *h_sum_hints, int *h_sum_lengths, int *h_sum_dirs, int no_sums)
{

    cout << "###h_sum_starts_x: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_starts_x[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_starts_y: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_starts_y[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_ends_x: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_ends_x[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_ends_y: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_ends_y[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_hints: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_hints[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_lengths: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_lengths[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_dirs: " << endl;
    for (int i = 0; i < no_sums; i++)
    {
        cout << h_sum_dirs[i] << " ";
    }
    cout << endl;
}

void flatten_sol_mat(int **sol_mat, int *h_sol_mat, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            h_sol_mat[i * n + j] = sol_mat[i][j];
        }
    }
}

void print_flattened_matrix(int *h_sol_mat, int m, int n)
{

    cout << "###Flattened matrix: " << endl;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << h_sol_mat[i * n + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

///////////////////////////////////////////////////
// Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////

///////////////////
// CUDA FUNCTIONS //
///////////////////

// For debugging.
__device__ void print_device_matrix(int **mat, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf(" %d", mat[i][j]);
        }
        printf("\n");
    }
}

// For debugging.
__device__ void print_device_matrix(int *mat, int size)
{
    for (int i = 0; i < size; i++)
    {
        printf(" %d", mat[i]);
    }
}

__device__ bool checkSumStatus(int remaining_sum, int remaining_cells)
{
    int current_max_num = 9;
    int current_min_num = 1;
    int max_num = 0;
    int min_num = 0;

    for (int i = 0; i < remaining_cells; i++)
    {
        max_num += current_max_num;
        min_num += current_min_num;
        current_max_num--;
        current_min_num++;
    }

    // remaining_sum > maximum value that can fit into remaining_cells.
    // We need to put bigger values to cells: anything containing smaller nums will be wrong
    if (remaining_sum > max_num)
        return false;

    // remaining_sum < minimum value that can fit into remaining_cells.
    // We need to put smaller values to cells: anything containing bigger nums will be wrong
    if (remaining_sum < min_num)
        return false;

    return true;
}

// Checks the solution matrix whether it is valid or has potential to be valid for a given sum object.
// It also checks for duplicates.
__device__ bool checkSum(int *&d_sum_starts_x, int *&d_sum_starts_y, int *&d_sum_ends_x, int *&d_sum_ends_y,
                         int *&d_sum_hints, int *&d_sum_dirs, int *&board, int m,
                         int n, int d_sum_idx, int k)
{
    // If there is no sum, don't check.
    if (d_sum_idx == -1)
    {
        return false;
    }

    int hint = d_sum_hints[d_sum_idx];

    int row_idx = d_sum_starts_x[d_sum_idx];
    int col_idx = d_sum_starts_y[d_sum_idx];

    // Check for a row sum
    if (d_sum_dirs[d_sum_idx] == 1)
    {
        int end_idx = d_sum_ends_y[d_sum_idx];

        // Continue iteration until there is a currently empty cell or end of the sum region.
        while (col_idx < end_idx && board[(row_idx * m) + col_idx] > 0)
        {
            // Substract the remaining sum by the value inside the sum region.
            hint -= board[(row_idx * m) + col_idx];
            bool status = checkSumStatus(hint, end_idx - col_idx - 1);
            // If sum status is not valid, return false.
            if (!status)
            {
                return false;
            }

            // Check for duplicates.
            if ((row_idx * m) + col_idx != k && board[k] == board[(row_idx * m) + col_idx])
            {
                return false;
            }

            col_idx++;
        }
    }

    // Check for a column sum
    else
    {
        int end_idx = d_sum_ends_x[d_sum_idx];

        // Continue iteration until there is a currently empty cell or end of the sum region.

        while (row_idx < end_idx && board[(row_idx * m) + col_idx] > 0)
        {
            // Substract the remaining sum by the value inside the sum region.
            hint -= board[(row_idx * m) + col_idx];
            bool status = checkSumStatus(hint, end_idx - row_idx - 1);

            // If sum status is not valid, return false.
            if (!status)
            {
                return false;
            }

            // Check for duplicates.
            if ((row_idx * m) + col_idx != k && board[k] == board[(row_idx * m) + col_idx])
            {

                return false;
            }

            row_idx++;
        }
    }
    return true;
}

// 2D array to map board cells to the flattened sum array indexes they are included in.
__device__ int **setCell2SumIdx(int *&d_sum_starts_x, int *&d_sum_starts_y, int *&d_sum_ends_x, int *&d_sum_ends_y,
                                int *&d_sum_dirs, int m, int n, int d_sum_count)
{
    int **cell_2_sum_idx = new int *[m * n];

    for (int i = 0; i < m * n; i++)
    {
        cell_2_sum_idx[i] = new int[2];
        for (int j = 0; j < 2; j++)
        {
            cell_2_sum_idx[i][j] = -1;
        }
    }

    for (int i = 0; i < d_sum_count; i++)
    {
        int start_row = d_sum_starts_x[i];
        int start_col = d_sum_starts_y[i];
        int end_row = d_sum_ends_x[i];
        int end_col = d_sum_ends_y[i];

        int start_k = start_row * m + start_col;
        int end_k = end_row * m + end_col;

        if (d_sum_dirs[i] == direction::d_right)
        {
            for (int j = start_k; j < end_k; j++)
            {
                if (cell_2_sum_idx[j][0] == -1) // If first sum
                {
                    cell_2_sum_idx[j][0] = i;
                }
                else
                {
                    cell_2_sum_idx[j][1] = i;
                }
            }
        }
        else
        {
            for (int j = start_k; j < end_k; j += m)
            {
                if (cell_2_sum_idx[j][0] == -1) // If first sum
                {
                    cell_2_sum_idx[j][0] = i;
                }
                else
                {
                    cell_2_sum_idx[j][1] = i;
                }
            }
        }
    }

    return cell_2_sum_idx;
}

// Generate deep copy of a matrix.
__device__ int **copyMatrix(int **mat, int m, int n)
{
    int **copy = new int *[m];
    for (int i = 0; i < m; i++)
    {
        copy[i] = new int[n];
        for (int j = 0; j < n; j++)
        {
            copy[i][j] = mat[i][j];
        }
    }

    return copy;
}

// Generate deep copy of a flatted matrix.
__device__ int *copyMatrixFlattened(int *mat, int size)
{
    int *copy = new int[size];
    for (int i = 0; i < size; i++)
    {
        copy[i] = mat[i];
    }
    return copy;
}

// Task generator kernel.
// Generates possibly valid boards from current tasks i.e. current boards.
// Each block is responsible for generating and checking 9 boards from previous boards. So, block_size = num_tasks.
// Each thread is responsible for a single board by inserting threadIdx.x + 1 to the next cell in the board.
// Store them in 2d array: new_tasks.
__global__ void kakuro_solver(int *d_sum_starts_x, int *d_sum_starts_y, int *d_sum_ends_x, int *d_sum_ends_y,
                              int *d_sum_hints, int *d_sum_dirs, int *d_sol_mat, int **tasks,
                              int m, int n, int k, int **d_cell2sum_idx, int **new_tasks, int dim)
{
    // Copy a board regarding to blockIdx.
    int *board = copyMatrixFlattened(tasks[blockIdx.x], m * n);

    // Insert a new number to board regarding to threadIdx.
    int num = threadIdx.x + 1;
    board[k] = num;

    // Get sum indexes from the map.
    int sum_idx_1 = d_cell2sum_idx[k][0];
    int sum_idx_2 = d_cell2sum_idx[k][1];
    bool status;

    // Check for the first sum.
    status = checkSum(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_dirs, board, m, n, sum_idx_1, k);

    if (!status)
    {
        new_tasks[blockDim.x * blockIdx.x + threadIdx.x] = nullptr;
        delete[] board;
        return;
    }

    // Check for the first sum.
    status = checkSum(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_dirs, board, m, n, sum_idx_2, k);
    if (!status)
    {
        new_tasks[blockDim.x * blockIdx.x + threadIdx.x] = nullptr;
        delete[] board;
        return;
    }

    // If no errors are present, generate a new task from the current state of the board.
    new_tasks[blockDim.x * blockIdx.x + threadIdx.x] = board;
}

// Kakuro solver controller kernel.
// This kernel runs with a single block and single thread.
// Since kakuro_solver is called multiple times for each step, it needs a controller inside the GPU.
// So that memory transfer does not cost as much as controlling it from the CPU.
__global__ void kakuro_kernel(int *d_sum_starts_x, int *d_sum_starts_y, int *d_sum_ends_x, int *d_sum_ends_y,
                              int *d_sum_hints, int *d_sum_dirs, int *d_sol_mat, int m, int n, int no_sums)
{
    int **cell_2_sum_idx = setCell2SumIdx(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y,
                                          d_sum_dirs, m, n, no_sums);

    // Start with a single board.
    int num_tasks = 1;
    int **tasks = new int *[num_tasks];

    tasks[0] = copyMatrixFlattened(d_sol_mat, m * n);

    for (int k = 0; k < m * n; k++)
    {
        if (tasks[0][k] != -2)
            continue;

        // Maximum size of tasks that will be generated is 9 * num_tasks, allocate memory for that.
        int num_new_tasks = 9 * num_tasks;
        int **new_tasks = new int *[num_new_tasks];

#ifdef DEBUG
        printf("NUM BLOCKS STEP %d: %d\n", k, num_tasks);
#endif

        // Run kakuro_solver for the current tasks / boards.
        kakuro_solver<<<num_tasks, 9>>>(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y,
                                        d_sum_hints, d_sum_dirs, d_sol_mat, tasks, m, n, k, cell_2_sum_idx, new_tasks, num_tasks);

        // Wait for kakuro_solver to finish.
        cudaDeviceSynchronize();

        // Organize and reset tasks and new tasks:

        // Delete previous tasks:
        for (int i = 0; i < num_tasks; i++)
        {
            delete[] tasks[i];
        }
        delete[] tasks;

        // Count new number of tasks.
        num_tasks = 0;
        for (int i = 0; i < num_new_tasks; i++)
        {
            if (new_tasks[i])
                num_tasks++;
        }

        // Allocate memory for new tasks.
        tasks = new int *[num_tasks];
        int task_idx = 0;

        for (int i = 0; i < num_new_tasks; i++)
        {
            if (new_tasks[i])
            {
                tasks[task_idx] = new_tasks[i];
                task_idx++;
            }
        }
        delete[] new_tasks;
    }

    for (int i = 0; i < m * n; i++)
    {
        d_sol_mat[i] = tasks[0][i];
    }
#ifdef DEBUG
    printf("SOL HERE\n");
    print_device_matrix(d_sol_mat, m * n);
    printf("\n\n");
#endif
}

///////////////////
// CUDA FUNCTIONS //
///////////////////

// Write solution to file.
void sol_mat_flattened_to_file(int **mat, int *d_sol_mat, int m, int n, string fname)
{
    ofstream to_write(fname);
    to_write << m << " " << n << "\n";

    int *h_sol_mat_f = new int[m * n];
    cudaMemcpy(h_sol_mat_f, d_sol_mat, m * n * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat[i][j] != -2)
                to_write << mat[i][j] << " ";
            else
                to_write << h_sol_mat_f[i * m + j] << " ";
        }
        to_write << "\n";
    }

    to_write.close();
}

int main(int argc, char **argv)
{

    std::string filename(argv[1]);
    std::ifstream file;
    file.open(filename.c_str());

    int m, n;
    file >> m;
    file >> n;

    int **mat;
    read_matrix(mat, file, m, n);
    print_one_matrix(mat, m, n);

    int **sol_mat;
    convert_sol(mat, sol_mat, m, n);
    // print_one_matrix(sol_mat, m, n);

    vector<sum> sums = get_sums(mat, m, n);

    // CUDA
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("==prop== Running on device: %d -- %s \n", 0, prop.name);
    printf("==prop== #of SM -- %d \n", prop.multiProcessorCount);
    printf("==prop== Max Threads Per Block: -- %d \n", prop.maxThreadsPerBlock);

    int grid_dim = 1;
    int block_dim = 1;

    int no_sums = sums.size();

    // Flattening sums and matrix
    int *h_sum_starts_x = new int[no_sums];
    int *h_sum_starts_y = new int[no_sums];
    int *h_sum_ends_x = new int[no_sums];
    int *h_sum_ends_y = new int[no_sums];
    int *h_sum_hints = new int[no_sums];
    int *h_sum_lengths = new int[no_sums];
    int *h_sum_dirs = new int[no_sums];

    flatten_sums(sums, h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

    print_flattened(h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

    int *h_sol_mat;
    h_sol_mat = new int[m * n];
    flatten_sol_mat(sol_mat, h_sol_mat, m, n);

    print_flattened_matrix(h_sol_mat, m, n);

    // Declare device pointers and copy data into device
    int *d_sum_starts_x, *d_sum_starts_y, *d_sum_ends_x, *d_sum_ends_y, *d_sum_hints, *d_sum_lengths, *d_sum_dirs, *d_sol_mat, *d_t_mats;

    cudaMalloc(&d_sum_starts_x, no_sums * sizeof(int));
    cudaMalloc(&d_sum_starts_y, no_sums * sizeof(int));
    cudaMalloc(&d_sum_ends_x, no_sums * sizeof(int));
    cudaMalloc(&d_sum_ends_y, no_sums * sizeof(int));
    cudaMalloc(&d_sum_hints, no_sums * sizeof(int));
    cudaMalloc(&d_sum_lengths, no_sums * sizeof(int));
    cudaMalloc(&d_sum_dirs, no_sums * sizeof(int));
    cudaMalloc(&d_sol_mat, (m * n) * sizeof(int));
    cudaMalloc(&d_t_mats, (m * n * grid_dim * block_dim) * sizeof(int)); // Allocating invidual matrix for each GPU thread
    // You may use this array if you will implement a thread-wise solution

    cudaMemcpy(d_sum_starts_x, h_sum_starts_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_starts_y, h_sum_starts_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_ends_x, h_sum_ends_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_ends_y, h_sum_ends_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_hints, h_sum_hints, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_lengths, h_sum_lengths, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_dirs, h_sum_dirs, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);

    // ALLOCATE 8GB
    size_t rsize = 1024ULL * 1024ULL * 1024ULL * 8ULL;
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, rsize);

    kakuro_kernel<<<grid_dim, block_dim>>>(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints,
                                           d_sum_dirs, d_sol_mat, m, n, no_sums);
    cudaDeviceSynchronize();
    // CUDA

    int *h_sol_mat_f = new int[m * n];
    cudaMemcpy(h_sol_mat_f, d_sol_mat, m * n * sizeof(int), cudaMemcpyDeviceToHost);
    print_flattened_matrix(h_sol_mat_f, m, n);

    string fname = argv[1];
    fname = fname.substr(0, fname.length() - 7) + "_solution.kakuro";

    sol_mat_flattened_to_file(mat, d_sol_mat, m, n, fname);
    // Similiar to sol_mat, use hints from mat and values from d_sol_mat

    for (int i = 0; i < n; i++)
    {
        delete mat[i];
        delete sol_mat[i];
    }

    delete mat;
    delete sol_mat;

    delete h_sum_starts_x;
    delete h_sum_starts_y;
    delete h_sum_ends_x;
    delete h_sum_ends_y;
    delete h_sum_hints;
    delete h_sum_lengths;
    delete h_sum_dirs;
    delete h_sol_mat;

    cudaFree(d_t_mats);
    cudaFree(d_sum_starts_x);
    cudaFree(d_sum_starts_y);
    cudaFree(d_sum_ends_x);
    cudaFree(d_sum_ends_y);
    cudaFree(d_sum_hints);
    cudaFree(d_sum_lengths);
    cudaFree(d_sum_dirs);
    cudaFree(d_sol_mat);

    return 0;
}
