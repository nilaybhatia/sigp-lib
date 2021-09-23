#include "Signal.hpp"

#include <assert.h>

#include <algorithm>
#include <complex>
#include <deque>
#include <iostream>
#include <vector>

Signal::Signal() {
    *this = Signal(0, std::vector<int>{});
}

Signal::Signal(int origin_index, const std::vector<int>& vals) {
    this->origin_index = origin_index;
    this->vals = vals;
}

Signal Signal::shift(int k) {
    // std::cout << *this << "k = " << k;
    std::deque<int> new_seq((this->vals).begin(), (this->vals).end());
    int new_origin_index;
    if (k > 0) {
        new_origin_index = std::max(this->origin_index - k, 0);
        for (int i = 1; i <= k - this->origin_index; i++) new_seq.push_front(0);
    } else {
        new_origin_index = this->origin_index - k;
        for (int i = 1; i <= new_origin_index - (int)this->vals.size() + 1; i++)
            new_seq.push_back(0);
    }
    return Signal(new_origin_index,
                  std::vector<int>(new_seq.begin(), new_seq.end()));
}

Signal Signal::reverse() {
    std::vector<int> reversed(this->vals.rbegin(), this->vals.rend());
    int new_origin_index = this->vals.size() - 1 - this->origin_index;
    return Signal(new_origin_index, reversed);
}

Signal Signal::reverse_and_shift(int k) {
    Signal reversed = this->reverse();
    Signal reverse_shifted = reversed.shift(-k);
    return reverse_shifted;
}

Signal Signal::time_scale(int c) {
    // TODO: Support any rational like a/b
    int new_origin_index = this->origin_index;
    std::deque<int> new_vals;
    new_vals.push_back(this->vals[new_origin_index]);

    // Shrinks
    for (int i = new_origin_index + 1; i < this->vals.size(); i++) {
        if ((i - new_origin_index) % c == 0) new_vals.push_back(this->vals[i]);
    }
    for (int i = new_origin_index - 1; i >= 0; i++) {
        if ((new_origin_index - i) % c == 0) new_vals.push_front(this->vals[i]);
    }
    return Signal(new_origin_index,
                  std::vector<int>(new_vals.begin(), new_vals.end()));
}

Signal Signal::multiply_scalar(int scalar) {
    std::vector<int> new_vals((this->vals).begin(), (this->vals).end());
    for (int i = 0; i < new_vals.size(); i++) {
        new_vals[i] *= scalar;
    }
    return Signal(this->origin_index, new_vals);
}

// int Signal::operator[] (int index){
//     if(index < 0 or index >= this->vals.size()) return 0;
//     else return this->vals[index];
// }

Signal operator*(const Signal& sig1, const Signal& sig2) {
    int new_orgin_index = std::min(sig1.origin_index, sig2.origin_index);
    std::deque<int> new_vals;
    int i = sig1.origin_index, j = sig2.origin_index;

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
    return Signal(new_orgin_index,
                  std::vector<int>{new_vals.begin(), new_vals.end()});
}

Signal operator+(const Signal& sig1, const Signal& sig2) {
    int new_orgin_index = std::max(sig1.origin_index, sig2.origin_index);
    std::deque<int> new_vals;
    int i = sig1.origin_index, j = sig2.origin_index;

    while (i >= 0 || j >= 0) {
        int a = (i < 0 ? 0 : sig1.vals[i]);
        int b = (j < 0 ? 0 : sig2.vals[j]);
        new_vals.push_front(a + b);
        i--;
        j--;
    }

    i = sig1.origin_index + 1, j = sig2.origin_index + 1;
    while (i < sig1.vals.size() || j < sig2.vals.size()) {
        int a = (i >= sig1.vals.size() ? 0 : sig1.vals[i]);
        int b = (j >= sig2.vals.size() ? 0 : sig2.vals[j]);
        new_vals.push_back(a + b);
        i++;
        j++;
    }
    return Signal(new_orgin_index,
                  std::vector<int>(new_vals.begin(), new_vals.end()));
}

std::ostream& operator<<(std::ostream& os, const Signal& seq) {
    os << "The signal is \n";
    for (int ele : seq.vals) os << ele << " ";
    os << "\t";
    os << "with origin as " << seq.origin_index << "\n";
    return os;
}

std::vector<std::vector<std::complex<double>>> Signal::multiply_matrices(
    const std::vector<std::vector<std::complex<double>>>& A,
    const std::vector<std::vector<int>>& B) {
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
                // std::cout << A[i][k] << " " << B[k][j] << "\n";
                std::complex<double> a = A[i][k];
                std::complex<double> b(B[k][j], 0);
                // why does b as a simple int not work?
                std::complex<double> prod(a * b);
                product[i][j] += prod;
            }
        }
    }
    return product;
}

std::vector<std::vector<int>> Signal::get_matrix(const std::vector<int>& v1,
                                                 const std::vector<int>& v2) {
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

Signal Signal::linear_convolution(const Signal& other) {
    int new_origin_index = this->origin_index + other.origin_index;
    int n = this->vals.size(), m = other.vals.size();

    std::vector<std::vector<int>> matrix = get_matrix(this->vals, other.vals);

    std::vector<int> result;
    for (int i = 0; i < n; i++) {
        int sum = 0;
        for (int row = i, col = 0; row >= 0; row--, col++) {
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    for (int j = 1; j < m; j++) {
        int sum = 0;
        for (int row = n - 1, col = j; col < m; row--, col++) {
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    Signal convoluted(new_origin_index, result);
    return convoluted;
}

Signal Signal::circular_convolution(const Signal& other) {
    int l = this->vals.size(), m = other.vals.size();
    int new_origin_index = 0;
    int n = std::max(l, m);

    std::vector<int> larger, smaller;
    if (l >= m) {
        larger.assign(this->vals.begin(), this->vals.end());
        smaller.assign(other.vals.begin(), other.vals.end());
    } else {
        smaller.assign(this->vals.begin(), this->vals.end());
        larger.assign(other.vals.begin(), other.vals.end());
    }
    // n x n matrix
    std::vector<std::vector<int>> matrix(n, std::vector<int>(n));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            matrix[i][j] = larger[(i - j + n) % n];
        }
    }
    for (int i = 1; i <= abs(l - m); i++) {
        smaller.push_back(0);
    }

    assert(smaller.size() == n);
    std::vector<int> result;
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

std::vector<int> Signal::get_correlation_vals(
    const std::vector<std::vector<int>>& matrix) {
    int n = matrix.size();
    std::vector<int> result;
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

Signal Signal::auto_correlate() {
    int n = this->vals.size();
    std::vector<std::vector<int>> matrix = get_matrix(this->vals, this->vals);
    std::vector<int> vals = get_correlation_vals(matrix);
    int new_origin_index = vals.size() / 2;
    Signal correlated(new_origin_index, vals);
    return correlated;
}

Signal Signal::cross_correlate(const Signal& other) {
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
    std::vector<std::vector<int>> matrix = get_matrix(smaller, larger);
    std::vector<int> result = get_correlation_vals(matrix);
    int new_origin_index = result.size() / 2;
    Signal correlated(new_origin_index, result);
    return correlated;
}

std::vector<std::complex<double>> Signal::DFT() {
    const double pi = 3.14159265358979323846;
    int N = this->vals.size();
    std::vector<std::vector<int>> vals_as_column_vector;
    for (int val : this->vals) {
        vals_as_column_vector.push_back({val});
    }
    std::vector<std::vector<std::complex<double>>> kernel(
        N, std::vector<std::complex<double>>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::complex<double> arg(0, -2 * pi * i * j / N);
            kernel[i][j] = exp(arg);
        }
    }
    /*for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            std::cout << kernel[i][j] << " ";
        }
        std::cout << "\n";
    }*/
    auto result = multiply_matrices(kernel, vals_as_column_vector);
    std::vector<std::complex<double>> result_vals;
    for (int i = 0; i < N; i++) result_vals.push_back(result[i][0]);
    return result_vals;
}

ComplexSignal::ComplexSignal(const std::vector<std::complex<double>>& vals) {
    this->vals = vals;
}
unsigned int ComplexSignal::bit_reverse(unsigned int num,
                                        unsigned int max_bits) {
    int ans = 0;
    for (int i = 0; i < max_bits; i++) {
        ans = (ans << 1);
        ans = (ans | (num & 1));
        num = (num >> 1);
    }
    return ans;
}
unsigned int ComplexSignal::lg(unsigned int num) {
    unsigned int ans = 0;
    while (num > 0) {
        ans++;
        num >>= 1;
    }
    return ans;
}
ComplexSignal ComplexSignal::FFT() {
    /*Ref:
     * https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Pseudocode*/
    const double pi = 3.14159265358979323846;
    int N = this->vals.size();
    if (N == 1) {
        return *this;
    } else {
        std::vector<std::complex<double>> temp;
        // bit reversal

        // max_bits = log2(N-1)
        unsigned int max_bits = lg(N - 1);

        // std::cout << "N=" << N << "max_bits=" << max_bits << "\n";
        for (int k = 0; k < N; k++) {
            unsigned int index = bit_reverse(k, max_bits);
            // std::cout << k << " " << index << "\n";
            temp.push_back(this->vals[index]);
        }

        // for(auto ele : temp) std::cout << ele << " "; std::cout << "\n";

        std::vector<std::complex<double>> even_vals(temp.begin(),
                                                    temp.begin() + N / 2);
        std::vector<std::complex<double>> odd_vals(temp.begin() + N / 2,
                                                   temp.end());
        ComplexSignal even_terms(even_vals);
        ComplexSignal odd_terms(odd_vals);
        ComplexSignal G = even_terms.FFT();
        ComplexSignal H = odd_terms.FFT();
        // std::cout << "N = " << N << "\n";
        // std::cout << "G=" << G << "H=" << H;
        std::vector<std::complex<double>> result;
        for (int k = 0; k < N; k++) {
            // std::cout << k << " ";
            std::complex<double> arg(0, -2 * pi * k / N);
            std::complex<double> W_N_k = exp(arg);
            result.push_back(G.vals[k % G.vals.size()] +
                             W_N_k * H.vals[k % H.vals.size()]);
        }
        return ComplexSignal(result);
    }
}
std::ostream& operator<<(std::ostream& os, const ComplexSignal& sig) {
    os << "The complex signal is \n";
    for (auto ele : sig.vals)
        os << std::noshowpos << std::real(ele) << std::showpos << " "
           << std::imag(ele) << "j" << std::noshowpos << "\t"
           << "(Amplitude: " << abs(ele) << ", Phase: " << arg(ele) << ")\n";
    return os;
}