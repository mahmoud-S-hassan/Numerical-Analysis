#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

const int MAX_N = 10;

void jacobiIteration(double A[MAX_N][MAX_N], double b[MAX_N], double x[MAX_N], int n, double Error)
{
    double newX[MAX_N];

    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += A[i][j] * x[j];
            }
        }
        newX[i] = (b[i] - sum) / A[i][i];
    }

    for (int i = 0; i < n; i++) {
        x[i] = newX[i];
    }
}

void gaussSeidelIteration(double A[MAX_N][MAX_N], double b[MAX_N], double x[MAX_N], int n, double Error)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += A[i][j] * x[j];
            }
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

bool hasConverged(double x[MAX_N], double prevX[MAX_N], int n, double Error)
{
    for (int i = 0; i < n; i++) {
        if (abs(x[i] - prevX[i]) >= Error) {
            return false;
        }
    }
    return true;
}

bool checkSolvability(double A[MAX_N][MAX_N], int n)
{
    for (int i = 0; i < n; i++) {
        double absA = abs(A[i][i]);
        double absB = 0.0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                absB += abs(A[i][j]);
            }
        }
        if (absA < absB) {
            return false;
        }
    }
    return true;
}

void generateTable(double A[MAX_N][MAX_N], double b[MAX_N], double xJacobi[MAX_N], double xSeidel[MAX_N], int n, double Error)
{
    int iterations = 0;

    cout << setprecision(4);
    // cout << "+----------+-----------+-----------+-----------+-----------+-----------+-----------+\n";

    cout << "+";
    for (int j = 0; j <= 2 * (n); j++) {
        cout << "----------+";
    }
    std::cout << std::endl;

    cout << "|" << left << setw(10) << "Iteration";
    for (int j = 0; j < n; j++) {
        cout << left << "|" << setw(10) << " Jacobi x" + to_string(j + 1) << "";
    }
    for (int j = 0; j < n; j++) {
        cout << left << "|" << setw(10) << " Seidel x" + to_string(j + 1) << "";
    }
    cout << "|\n";

    // cout << "+----------+-----------+-----------+-----------+-----------+-----------+-----------+\n";
    cout << "+";
    for (int j = 0; j <= 2 * (n); j++) {
        cout << "----------+";
    }

    std::cout << std::endl;
    while (true) {
        cout << "|  " << setw(8) << iterations;

        double prevXJacobi[MAX_N];
        double prevXSeidel[MAX_N];

        for (int i = 0; i < n; i++) {
            prevXJacobi[i] = xJacobi[i];
            prevXSeidel[i] = xSeidel[i];
        }

        // Jacobi iteration
        jacobiIteration(A, b, xJacobi, n, Error);
        for (int j = 0; j < n; j++) {
            cout << "|" << setw(10) << xJacobi[j];
        }

        // Gauss-Seidel iteration
        gaussSeidelIteration(A, b, xSeidel, n, Error);
        for (int j = 0; j < n; j++) {
            cout << "|" << setw(10) << xSeidel[j];
        }
        cout << "|\n";
        // cout << "+----------+-----------+-----------+-----------+-----------+-----------+-----------+\n";

        cout << "+";
        for (int j = 0; j <= 2 * (n); j++) {
            cout << "----------+";
        }
        std::cout << std::endl;
        bool hasJacobiConverged = hasConverged(xJacobi, prevXJacobi, n, Error);
        bool hasSeidelConverged = hasConverged(xSeidel, prevXSeidel, n, Error);

        if (hasJacobiConverged && hasSeidelConverged) {
            break;
        }

        iterations++;

        // break loop #for testing

        // if (iterations == 70) {
        //     return;
        // }
    }
}

void rearrangeAndSolve(double A[MAX_N][MAX_N], double B[MAX_N], double xJacobi[MAX_N], double xSeidel[MAX_N], int n, double Error)
{
    int equationOrder[MAX_N];
    for (int i = 0; i < n; i++) {
        equationOrder[i] = i;
    }
    bool JacobiSeidel = false;
    double rearrangedA[MAX_N][MAX_N];
    double rearrangedB[MAX_N];

    do {
        // rearranging

        for (int i = 0; i < n; i++) {
            int rearrangedIndex = equationOrder[i];
            for (int j = 0; j < n; j++) {
                rearrangedA[i][j] = A[rearrangedIndex][j];
            }
            rearrangedB[i] = B[rearrangedIndex];
        }

        if (checkSolvability(rearrangedA, n)) {
            generateTable(rearrangedA, rearrangedB, xJacobi, xSeidel, n, Error);
            JacobiSeidel = true;
            return;
        }

    } while (next_permutation(equationOrder, equationOrder + n));

    if (!JacobiSeidel) {

        generateTable(A, B, xJacobi, xSeidel, n, Error);
    }
}

int main()
{
    double Error = 0.01;
    int N;
    cout << "Enter the matrix dimension (n): ";
    cin >> N;

    double A[MAX_N][MAX_N];
    double b[MAX_N];

    double xJacobi[MAX_N] = { 0.0 };
    double xSeidel[MAX_N] = { 0.0 };

    cout << "Enter the matrix elements:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << "A[" << i << "][" << j << "]: ";
            cin >> A[i][j];
        }
    }

    cout << "Enter the RHS vector:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "b[" << i << "]: ";
        cin >> b[i];
    }

    // --------------------------------------------
    // int N = 2;
    // double A[MAX_N][MAX_N] = {
    //     { 2, -3, 0 },
    //     { 1, 3, -10 },
    //     { 3, 0, 1 }
    // };
    //
    // double b[MAX_N] = { -7, 9, 13 };

    rearrangeAndSolve(A, b, xJacobi, xSeidel, N, Error);

    return 0;
}
