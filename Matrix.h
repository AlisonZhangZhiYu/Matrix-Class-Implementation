/*
 * default-initialized method -> separate the allocation
 * of memory and the construction of objects.
 * */

#ifndef MATRIX_H
#define MATRIX_H

#include <memory>
#include <iostream>
#include <string>
#include <random>
#include <type_traits>


template <typename T>
class Matrix;

template <typename T>
std::ostream &print_matrix(std::ostream &os, const Matrix<T> &matrix);

template <typename T>
class Matrix {
    friend std::ostream &print_matrix<>(std::ostream &os, const Matrix<T> &matrix);
public:
    Matrix();
    Matrix(size_t row, size_t col);
    Matrix(size_t row, size_t col, const T* input);
    ~Matrix();

    inline int get_row() const;
    inline int get_col() const;
    inline T get(size_t row, size_t col) const;
    inline T *get_ptr(size_t row, size_t col) const;
    inline void set(size_t row, size_t col, T target);
    inline size_t get_count() const;
    inline std::pair<size_t, size_t> get_dim();
    inline bool empty() const;

    // Calculation of matrix
    Matrix<T> operator+(const Matrix<T> &m1) const;
    Matrix<T> operator-(const Matrix<T> &m1) const;
    Matrix<T> operator*(const Matrix<T> &m1) const;
    Matrix<T> operator/(const Matrix<T> &m1) const;

    // operation with scalar
    Matrix<T> operator*(const T &s1) const;
    Matrix<T> operator/(const T &s1) const;

    // Inner Product
    T inner(const Matrix<T> &matrix) const;
    // Dot Product
    Matrix<T> dot(const Matrix<T> &matrix);
    // Transpose
    Matrix<T> transpose();

    Matrix<T> flatten();

    void random(const T &mean, const T &stdv);

    // Set to the identity Matrix
    void identity(const size_t &n);

    void resize(const size_t &row, const size_t &col);

    void zeros(const size_t &n);

    void zeros(const size_t &row, const size_t &col);

    void ones(const size_t &n);

    void ones(const size_t &row, const size_t &col);


private:
    std::pair<size_t, size_t> dim; // row, col
    size_t count{};
    T* data_ptr;

    inline std::string get_shape() const;
};

template<typename T>
std::string Matrix<T>::get_shape() const {
    return std::string("(" + std::to_string(this->dim.first) + "," + std::to_string(this->dim.second) + ")");
}

/*
 * Type: int
 * [[1, 2, 3],
 *  [4, 5, 6],
 *  [7, 8, 9]]
 *  [1 2 3 4]
 *
 * */
template <typename T>
std::ostream &print_matrix(std::ostream &os, const Matrix<T> &matrix) {
    // if there's no element in the matrix, then return 0;
    if (matrix.count == 0) {
        os << 0;
        return os;
    }

    // One dimension or Two dimension?
    bool flag = false;
    if (matrix.dim.first != 0 && matrix.dim.second != 0) {
        flag = true;
        os << "[";
    }

    os << "[";
    for (size_t i = 0; i < matrix.count; i++) {
        os << *(matrix.data_ptr + i);
        if ((i + 1) % matrix.get_col() != 0)
            os << " ";
        else {
            os << "]";
            if (i + 1 != matrix.count)
                os << std::endl << " [";
        }
    }

    if (flag)
        os << "]" << std::endl;

    return os;
}

template<typename T>
Matrix<T>::~Matrix() {
    delete[] data_ptr;
}

template<typename T>
Matrix<T>::Matrix() {
    this->dim.first = 0;
    this->dim.second = 0;
    this->count = 0;
}

template<typename T>
Matrix<T>::Matrix(size_t row, size_t col) {
    this->dim.first = row;
    this->dim.second = col;
    count = (row == 0 || col == 0) ? (row == 0 ? col : row) : row * col;

    // default initialize the matrix
    data_ptr = new T[count];
    for (T *p = data_ptr; p != data_ptr + count; p++)
        *p = 0.0;
}

template<typename T>
Matrix<T>::Matrix(size_t row, size_t col, const T *input) {
    this->dim.first = row;
    this->dim.second = col;
    count = (row == 0 || col == 0) ? (row == 0 ? col : row) : row * col;

    data_ptr = new T[count];
    for (T *p = data_ptr; p != data_ptr + count; p++)
        *p = *input++;
}

template<typename T>
int Matrix<T>::get_col() const {
    return dim.second;
}

template<typename T>
int Matrix<T>::get_row() const {
    return dim.first;
}

template<typename T>
T Matrix<T>::get(size_t row, size_t col) const {
    if (data_ptr == nullptr || empty())
        throw std::runtime_error("Nullptr for matrix! Matrix hasn't been initialized! \n Initialize first before using it!");

    return *(data_ptr + (row * this->dim.second + col));
}

template<typename T>
T *Matrix<T>::get_ptr(size_t row, size_t col) const {
    if (data_ptr == nullptr || empty())
        throw std::runtime_error("Nullptr for matrix! Matrix hasn't been initialized! \n Initialize first before using it!");

    return (data_ptr + (row * this->dim.second + col));
}

template<typename T>
void Matrix<T>::set(size_t row, size_t col, T target) {
    if (data_ptr == nullptr || empty())
        throw std::runtime_error("Nullptr for matrix! Matrix hasn't been initialized! \n Initialize first before using it!");

    if (row >= dim.first || col >= dim.second)
        throw std::out_of_range("out of dimension error!");
    *(data_ptr + (row * dim.second + col)) = target;
}

template<typename T>
std::pair<size_t, size_t> Matrix<T>::get_dim() {
    return this->dim;
}

template<typename T>
size_t Matrix<T>::get_count() const {
    return this->count;
}

template<typename T>
bool Matrix<T>::empty() const {
    return this->count == 0;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &m1) const {
    if (this->dim.first != m1.dim.first || this->dim.second != m1.dim.second)
        throw std::invalid_argument("Matrix dimension must agree!");

    Matrix<T> result(this->dim.first, this->dim.second);
    for (size_t i = 0; i < this->count; i++) {
        *(result.data_ptr + i) = *(this->data_ptr + i) + *(m1.data_ptr + i);
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &m1) const {
    if (this->dim.first != m1.dim.first || this->dim.second != m1.dim.second)
        throw std::invalid_argument("Matrix dimension must agree!");

    Matrix<T> result(this->dim.first, this->dim.second);
    for (size_t i = 0; i < this->count; i++) {
        *(result.data_ptr + i) = *(this->data_ptr + i) - *(m1.data_ptr + i);
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &m1) const {
    if (this->dim.first != m1.dim.first || this->dim.second != m1.dim.second)
        throw std::invalid_argument("Matrix dimension must agree!");

    Matrix<T> result(this->dim.first, this->dim.second);
    for (size_t i = 0; i < this->count; i++) {
        *(result.data_ptr + i) = *(this->data_ptr + i) * *(m1.data_ptr + i);
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const Matrix<T> &m1) const {
    if (this->dim.first != m1.dim.first || this->dim.second != m1.dim.second)
        throw std::invalid_argument("Matrix dimension must agree!");

    Matrix<T> result(this->dim.first, this->dim.second);
    for (size_t i = 0; i < this->count; i++) {
        if (*(m1.data_ptr + i) == 0)
            throw std::runtime_error("Runtime Warning: divided by zero!")
        *(result.data_ptr + i) = *(this->data_ptr + i) / *(m1.data_ptr + i);
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const T &s1) const {
    if (s1 == 0)
        throw std::runtime_error("RuntimeWarning: divided by scalar zero!");

    Matrix<T> result(this->dim.first, this->dim.second);
    for (size_t i = 0; i < this->count; i++) {
        *(result.data_ptr + i) = *(this->data_ptr + i) / s1;
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T &s1) const {
    Matrix<T> result(this->dim.first, this->dim.second);
    for (size_t i = 0; i < this->count; i++) {
        *(result.data_ptr + i) = *(this->data_ptr + i) * s1;
    }

    return result;
}

template<typename T>
T Matrix<T>::inner(const Matrix<T> &m1) const {
    if (this->dim.first != m1.dim.first || this->dim.second != m1.dim.second)
        throw std::invalid_argument("Matrix dimension must agree!");

    T res = 0;
    for (size_t i = 0; i < count; i++)
        res += *(m1.data_ptr + i) * *(this->data_ptr + i);

    return res;
}

/*
 * 当前为按照行存储
 * 按列存储理论上更高效
 *
 * */

template<typename T>
Matrix<T> Matrix<T>::dot(const Matrix<T> &mat) {
    if (this->dim.second != mat.dim.first) {
        std::string msg("Invalid Argument: Shape " + get_shape() + " and " + mat.get_shape() + " are not aligned!");
        throw std::invalid_argument(msg.c_str());
    }

    Matrix<T> res(this->get_row(), mat.get_col());

    int idx = 0;
    for (size_t i = 0; i < this->get_row(); i++) {
        for (size_t j = 0; j < mat.get_col(); j++) {
            for (size_t k = 0; k < this->get_col(); k++)
                *(res.data_ptr + idx) += this->get(i, k) * mat.get(k, j);
            idx++;
        }
    }

    return res;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() {
    if (data_ptr == nullptr || empty())
        throw std::runtime_error("Nullptr for matrix! Matrix hasn't been initialized! \n Initialize first before using it!");

    Matrix<T> mat(this->dim.second, this->dim.first);
    for (size_t i = 0; i < this->dim.second; i++) {     // col
        for (size_t j = 0; j < this->dim.first; j++) {  // row
            *(mat.data_ptr + i * this->dim.first + j) = this->get(j, i);
        }
    }

    return mat;
}

template<typename T>
Matrix<T> Matrix<T>::flatten() {
    Matrix<T> res(0, this->count, this->data_ptr);
    return res;
}

// stdv: standard deviation
template<typename T>
void Matrix<T>::random(const T &mean, const T &stdv) {
    std::random_device rd; // 随机数生成器种子
    std::mt19937 gen(rd()); // Mersenne Twister 引擎

    if constexpr (std::is_floating_point<T>::value) {
        std::normal_distribution<T> dist(mean, stdv); // 正态分布
        for (size_t i = 0; i < count; i++)
            *(data_ptr + i) = dist(gen);

    } else if constexpr (std::is_integral<T>::value) {
        std::uniform_int_distribution<T> dist(mean, stdv); // 均匀分布
        for (size_t i = 0; i < count; i++)
            *(data_ptr + i) = dist(gen);

    } else {
        static_assert(std::is_arithmetic<T>::value, "Template argument must be an arithmetic type.");
    }
}

template<typename T>
void Matrix<T>::identity(const size_t &n) {
    if (dim.first != n || dim.second != n) {
        // resize the matrix
        resize(n, n);
    }

    // update the matrix data
    int delta = 0;
    for (T *ptr = data_ptr; ptr != data_ptr + count; ptr+=n) {
        *(ptr + delta) = 1;
        delta++;
    }
}

template<typename T>
void Matrix<T>::resize(const size_t &row, const size_t &col) {
    if (dim.first == row && dim.second == col)
        return;

    if (!empty())
        delete[] data_ptr;

    // Update the dimension of matrix
    this->dim.first = row;
    this->dim.second = col;
    this->count = (row == 0 || col == 0) ? (row == 0 ? col : row) : row * col;

    // Reallocate Memory (Improve by using allocator)
    data_ptr = new T[this->count];

    for (T *ptr = data_ptr; ptr != data_ptr + count; ptr++)
        *ptr = 0.0;
}

template<typename T>
void Matrix<T>::zeros(const size_t &row, const size_t &col) {
    if (empty() || dim.first != row || dim.second != col)
        resize(row, col);

    for (T *ptr = data_ptr; ptr != data_ptr + count; ptr++)
        *ptr = 0.0;
}

template<typename T>
void Matrix<T>::zeros(const size_t &n) {
    zeros(n, n);
}

template<typename T>
void Matrix<T>::ones(const size_t &row, const size_t &col) {
    if (empty() || dim.first != row || dim.second != col)
        resize(row, col);

    for (T *ptr = data_ptr; ptr != data_ptr + count; ptr++)
        *ptr = 1.0;
}

template<typename T>
void Matrix<T>::ones(const size_t &n) {
    ones(n, n);
}


#endif
