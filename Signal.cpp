#include "Signal.hpp"

#include <assert.h>

#include <algorithm>
#include <complex>
#include <deque>
#include <iostream>
#include <stdexcept>
#include <vector>

template <class T>
Signal<T>::Signal() {
    *this = Signal(0, std::vector<T>{});
}

template <class T>
Signal<T>::Signal(T origin_index, const std::vector<T>& vals) {
    this->origin_index = origin_index;
    this->vals = vals;
}

template <class T>
Signal<T> Signal<T>::shift(int k) {
    // std::cout << *this << "k = " << k;
    std::deque<T> new_seq((this->vals).begin(), (this->vals).end());
    T new_origin_index;
    if (k > 0) {
        new_origin_index = std::max(this->origin_index - k, 0);
        for (int i = 1; i <= k - this->origin_index; i++) new_seq.push_front(0);
    } else {
        new_origin_index = this->origin_index - k;
        for (int i = 1; i <= new_origin_index - (int)this->vals.size() + 1; i++)
            new_seq.push_back(0);
    }
    return Signal<T>(new_origin_index,
                     std::vector<T>(new_seq.begin(), new_seq.end()));
}

template <class T>
Signal<T> Signal<T>::reverse() {
    std::vector<T> reversed(this->vals.rbegin(), this->vals.rend());
    T new_origin_index = this->vals.size() - 1 - this->origin_index;
    return Signal<T>(new_origin_index, reversed);
}

template <class T>
Signal<T> Signal<T>::reverse_and_shift(int k) {
    Signal reversed = this->reverse();
    Signal reverse_shifted = reversed.shift(-k);
    return reverse_shifted;
}

template <class T>
Signal<T> Signal<T>::time_scale(int c) {
    // TODO: Support any rational like a/b
    int new_origin_index = this->origin_index;
    std::deque<T> new_vals;
    new_vals.push_back(this->vals[new_origin_index]);

    // Shrinks
    for (int i = new_origin_index + 1; i < this->vals.size(); i++) {
        if ((i - new_origin_index) % c == 0) new_vals.push_back(this->vals[i]);
    }
    for (int i = new_origin_index - 1; i >= 0; i++) {
        if ((new_origin_index - i) % c == 0) new_vals.push_front(this->vals[i]);
    }
    return Signal(new_origin_index,
                  std::vector<T>(new_vals.begin(), new_vals.end()));
}

template <class T>
Signal<T> Signal<T>::multiply_scalar(int scalar) {
    std::vector<T> new_vals((this->vals).begin(), (this->vals).end());
    for (int i = 0; i < new_vals.size(); i++) {
        new_vals[i] *= scalar;
    }
    return Signal(this->origin_index, new_vals);
}

// int Signal::operator[] (int index){
//     if(index < 0 or index >= this->vals.size()) return 0;
//     else return this->vals[index];
// }
template <class T1, class T2>
auto operator*(const Signal<T1>& sig1, const Signal<T2>& sig2) {
    T new_orgin_index = std::min(sig1.origin_index, sig2.origin_index);
    std::deque<T> new_vals;
    T i = sig1.origin_index, j = sig2.origin_index;

    while (i >= 0 && j >= 0) {
        new_vals.push_front(sig1.vals[i] * sig2.vals[j]);
        i--;
        j--;
    }

    i = sig1.origin_index + 1, j = sig2.origin_index + 1;
    while (i < sig1.vals.size() && j < sig2.vals.size()) {
        new_vals.push_back(sig1.vals[i] * sig2.vals[j]);
        i++;
        j++;
    }
    return Signal<T>(new_orgin_index,
                     std::vector<T>{new_vals.begin(), new_vals.end()});
}

template <class T>
Signal<T> operator+(const Signal<T>& sig1, const Signal<T>& sig2) {
    T new_orgin_index = std::max(sig1.origin_index, sig2.origin_index);
    std::deque<T> new_vals;
    T i = sig1.origin_index, j = sig2.origin_index;

    while (i >= 0 || j >= 0) {
        T a = (i < 0 ? 0 : sig1.vals[i]);
        T b = (j < 0 ? 0 : sig2.vals[j]);
        new_vals.push_front(a + b);
        i--;
        j--;
    }

    i = sig1.origin_index + 1, j = sig2.origin_index + 1;
    while (i < sig1.vals.size() || j < sig2.vals.size()) {
        T a = (i >= sig1.vals.size() ? 0 : sig1.vals[i]);
        T b = (j >= sig2.vals.size() ? 0 : sig2.vals[j]);
        new_vals.push_back(a + b);
        i++;
        j++;
    }
    return Signal<T>(new_orgin_index,
                     std::vector<T>(new_vals.begin(), new_vals.end()));
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Signal<T>& sig) {
    os << "The signal is \n";
    // if (typeid(T).name() == std::complex) {
    //     for (T ele : sig.vals)
    //         os << std::real(ele) << " " << std::imag(ele) << "j";
    //     os << "\t";
    //     os << "with origin as " << sig.origin_index << "\n";
    // } else {
        for (T ele : sig.vals) os << ele << " ";
        os << "\t";
        os << "with origin as " << sig.origin_index << "\n";
    //}
    return os;
}

// very makeshift, should probably use eigen sth
template <class T>
std::vector<std::vector<T>> Signal<T>::multiply_matrices(
    const std::vector<std::vector<T>>& A,
    const std::vector<std::vector<T>>& B) {
    if (A[0].size() != B.size()) {
        throw std::invalid_argument("Matrices are not compatible");
    }
    int m = A.size();
    int n = A[0].size();
    int p = B[0].size();

    std::vector<std::vector<std::complex<double>>> product(
        m, std::vector<std::complex<double>>(p));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                product[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return product;
}
template <class T>
std::vector<std::vector<T>> Signal<T>::get_matrix(const std::vector<T>& v1,
                                                  const std::vector<T>& v2) {
    // returns the shortcut product matrix i.e the matrix
    // formed by tabular product of the 2 signals
    // v1 vertical and v2 horizontal
    int n = v1.size(), m = v2.size();

    // n x m matrix
    std::vector<std::vector<int>> matrix(n, std::vector<int>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = v1[i] * v2[j];
        }
    }
    return matrix;
}

template <class T>
Signal<T> Signal<T>::linear_convolution(const Signal<T>& other) {
    T new_origin_index = this->origin_index + other.origin_index;
    int n = this->vals.size(), m = other.vals.size();

    std::vector<std::vector<T>> matrix = get_matrix(this->vals, other.vals);

    std::vector<T> result;
    for (int i = 0; i < n; i++) {
        T sum = 0;
        for (int row = i, col = 0; row >= 0; row--, col++) {
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    for (int j = 1; j < m; j++) {
        T sum = 0;
        for (int row = n - 1, col = j; col < m; row--, col++) {
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    Signal convoluted(new_origin_index, result);
    return convoluted;
}

template <class T>
Signal<T> Signal<T>::circular_convolution(const Signal<T>& other) {
    int l = this->vals.size(), m = other.vals.size();
    T new_origin_index = 0;
    int n = std::max(l, m);

    std::vector<T> larger, smaller;
    if (l >= m) {
        // TODO: maybe use const references?
        larger.assign(this->vals.begin(), this->vals.end());
        smaller.assign(other.vals.begin(), other.vals.end());
    } else {
        smaller.assign(this->vals.begin(), this->vals.end());
        larger.assign(other.vals.begin(), other.vals.end());
    }
    // n x n matrix
    std::vector<std::vector<T>> matrix(n, std::vector<T>(n));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            matrix[i][j] = larger[(i - j + n) % n];
        }
    }
    for (int i = 1; i <= abs(l - m); i++) {
        smaller.push_back(0);
    }

    assert(smaller.size() == n);
    std::vector<T> result;
    for (int row = 0; row < n; row++) {
        int sum = 0;
        for (int col = 0; col < n; col++) {
            sum += matrix[row][col] * smaller[col];
        }
        result.push_back(sum);
    }

    Signal convoluted(new_origin_index, result);
    return convoluted;
}

template <class T>
std::vector<T> Signal<T>::get_correlation_vals(
    const std::vector<std::vector<T>>& matrix) {
    int n = matrix.size();
    std::vector<T> result;
    for (int i = n - 1; i >= 0; i--) {
        int sum = 0;
        for (int row = i, col = 0; row < n; row++, col++) {
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    for (int j = 1; j < n; j++) {
        int sum = 0;
        for (int row = 0, col = j; col < n; row++, col++) {
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    return result;
}

template <class T>
Signal<T> Signal<T>::auto_correlate() {
    int n = this->vals.size();
    std::vector<std::vector<T>> matrix = get_matrix(this->vals, this->vals);
    std::vector<T> vals = get_correlation_vals(matrix);
    int new_origin_index = vals.size() / 2;
    Signal correlated(new_origin_index, vals);
    return correlated;
}

template <class T>
Signal<T> Signal<T>::cross_correlate(const Signal& other) {
    int l = this->vals.size(), m = other.vals.size();
    int n = std::max(l, m);

    std::vector<int> larger, smaller;
    if (l >= m) {
        larger.assign(this->vals.begin(), this->vals.end());
        smaller.assign(other.vals.begin(), other.vals.end());
    } else {
        smaller.assign(this->vals.begin(), this->vals.end());
        larger.assign(other.vals.begin(), other.vals.end());
    }
    for (int i = 1; i <= abs(l - m); i++) {
        smaller.push_back(0);
    }

    assert(smaller.size() == n);
    // n x n matrix
    // TODO: Since cross correlation is not commutative, order is important
    // Rn, smaller signal must be the 2nd one
    std::vector<std::vector<T>> matrix = get_matrix(smaller, larger);
    std::vector<T> result = get_correlation_vals(matrix);
    T new_origin_index = result.size() / 2;
    Signal correlated(new_origin_index, result);
    return correlated;
}

template <class T>
Signal<T> Signal<T>::DFT() {
    const double pi = 3.14159265358979323846;
    int N = this->vals.size();
    std::vector<std::vector<T>> vals_as_column_vector;
    for (T val : this->vals) {
        vals_as_column_vector.push_back({val});
    }
    std::vector<std::vector<std::complex<double>>> kernel(
        N, std::vector<std::complex<double>>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            kernel[i][j] = exp(-2 * pi * i * j / N);
        }
    }
    auto result = multiply_matrices(kernel, vals_as_column_vector);
    std::vector<std::complex<double>> result_vals;
    for (int i = 0; i < N; i++) result_vals.push_back(result[i][0]);
    return Signal(0, result_vals);
}