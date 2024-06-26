#include <iostream>
#include <string>
#include <chrono>
#include "Matrix.h"

using namespace std;

int main() {
    int *arr = new int[100]();
    for (int i = 0; i < 100; i++) {
        arr[i] = i;
    }

    Matrix<int> mat1(10, 10, arr);
    Matrix<int> mat2(10, 10, arr);
    Matrix<int> res = mat1.dot(mat2);

}
