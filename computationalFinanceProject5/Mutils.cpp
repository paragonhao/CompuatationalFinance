//
// Created by paragonhao on 16/1/19.
//

#include "Mutils.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

double Mutils::Mean(double arr[], int n)
{
    double mean = 0.0;
    for (int idx = 0; idx < n; idx++)
    {
        mean += arr[idx];
    }
    mean /= static_cast<double>(n);
    return mean;
}


double Mutils::StDev(double arr[], int n)
{
    double mean = Mutils::Mean(arr, n);
    double variance = 0.0;
    for (size_t idx = 0; idx < n; idx++)
    {
        double temp = arr[idx] - mean;
        variance += temp*temp;
    }

    // Compute sample variance using Bessel's correction (see http://en.wikipedia.org/wiki/Bessel%27s_correction)
    variance /= static_cast<double>(n) - (n == 1 ? 0.0 : 1.0);

    // Standard deviation is square root of variance
    return std::sqrt(variance);
}

double Mutils::Corr(double *x, double *y, int size)
{
    double meanX = Mutils::Mean(x, size);
    double meanY = Mutils::Mean(y, size);
    double numerator = 0.0;
    double denominatorX = 0.0;
    double denominatorY = 0.0;

    for(int i = 0; i< size; i++){
        numerator += (x[i] - meanX) * (y[i] - meanY);
        denominatorX += (x[i] - meanX) * (x[i] - meanX);
        denominatorY += (y[i] - meanY) * (y[i] - meanY);
    }

    numerator /= (size-1);
    denominatorX = sqrt(denominatorX /(size - 1));
    denominatorY = sqrt(denominatorY /(size -1));

    return numerator/(denominatorY* denominatorX);
}

double Mutils::Cov(double *x, double *y, int size){
    double meanX = Mutils::Mean(x, size);
    double meanY = Mutils::Mean(y, size);
    double sse = 0 ;
    for(int i = 0; i< size; i++){
        sse += (x[i] - meanX) * (y[i] - meanY);
    }
    return sse/(size-1);
}



void Mutils::WriteArrayToCSV( double *arr, int n, const string filename){
    ofstream file;
    file.open(filename);
    for(int i=0; i < n; i++){
        file << arr[i] << ",\n";
    }
    file.close();
}


void Mutils::WriteToCSV2DMatrix(double arr[][31], int row,int col, const string filename){
    ofstream file;
    file.open(filename);

    for(int i=0; i<row;i++){
        for(int j=0; j<col; j++){
            file << arr[j][i] << ",";
        }
        file << "\n";
    }

    file.close();
}

double * Mutils::MatrixMultiply(double *arr, int size, double t){

    for(int i=0; i< size; i++){
        arr[i] *= t;
    }

    return arr;
}

double * Mutils::MatrixAddition(double *arr, int size, double num){

    for(int i=0; i< size; i++){
        arr[i] += num;
    }

    return arr;
}

double Mutils::pnorm(double x)
{
    const double d1 = 0.0498673470;
    const double d2 = 0.0211410061;
    const double d3 = 0.0032776263;
    const double d4 = 0.0000380036;
    const double d5 = 0.0000488906;
    const double d6 = 0.0000053830;

    double probability = 0.0;
    double absX = (x >=0)?x:-x;
    double temp =  1 + (d1 * absX) + (d2 * absX * absX) + (d3 * absX * absX * absX) + (d4 * pow(absX,4))
            + (d5 * pow(absX,5)) + (d6 * pow(absX, 6));
    probability = 1 - 0.5 * pow(temp, -16);

    return  (x >= 0)? probability : (1- probability);
}


/* Vector element wise multiplication
 *
 * @param:
 * double *: first array
 * double *: seconed array
 * int : size of the array to check
 * */
double Mutils::arrayElementWiseMultiply(double * v1, double * v2, int NumBits){
    double sum = 0;
    for(int i = 0; i<NumBits; i++){
        sum += v1[i] * v2[i];
    }
    return sum;
}

/* Get cofator
 *
 * @param:
 * vector<vector<double>> : original matrix
 * vector<vector<double>> : temp
 * int p, q: row and col to get a cofactor
 * int n: n is current dimension of vector matrix
 * */
void getCofactor(vector<vector<double>> &matrix, vector<vector<double>> &temp, int p, int q, int n)
{
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = matrix[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
 *
 * @param:
 * double ** A: matrix
 * int order: order of matrix
 *
*/
double determinant(vector<vector<double>> &matrix, int n)
{
    int D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return matrix[0][0];

    vector<vector< double> > temp (n, vector<double>(n, 0)); // To store cofactors

    int sign = 1;  // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(matrix, temp, 0, f, n);
        D += sign * matrix[0][f] * determinant(temp, n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

/* Function to get adjoint of A[N][N] in adj[N][N].
 *
 * @param:
 * double ** A: original matrix
 * double ** adj: adjoin matrix
 * int order: order of matrix
 *
*/
void adjoint(vector<vector<double>> &matrix, vector<vector<double>> &adj, int N)
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;
//    vector<double> temp;
//    double * temp[N];

    vector<vector< double> > temp (N, vector<double>(N, 0));

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(matrix, temp, i, j, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinant(temp, N-1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(vector<vector<double>> &matrix, vector<vector<double>> &inverse, int N)
{
    // Find determinant of A[][]
    double det = determinant(matrix, N);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    vector<vector< double> > adj (N, vector<double>(N, 0));
    adjoint(matrix, adj, N);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i][j] = adj[i][j]/det;

    return true;
}

/* Matrix Inversion
 * A simple program adapted from: https://www.geeksforgeeks.org/adjoint-inverse-matrix/
 *
 * @param
 *  vector< vector< double > >: pass in a matrix
 *
 * */
vector<vector<double>> Mutils::MatrixInverse(vector<vector<double>> &matrix, int order) {
    vector< vector< double > > adj ( order, vector<double> ( order, 0 ) );
    vector< vector< double > > invMatrix ( order, vector<double> ( order, 0 ) );

    adjoint(matrix, adj, order);
    inverse(matrix, invMatrix,order);

    return invMatrix;
}


double Mutils::max(const double & num1, const double & num2){
    if(num1 >= num2){
        return num1;
    }else{
        return num2;
    }

}


