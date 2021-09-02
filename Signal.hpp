/*
Author: Nilay Bhatia (@nilaybhatia)
Date: 13th August 2021
References:
    https://docs.microsoft.com/en-us/cpp/cpp/header-files-cpp?view=msvc-160
    https://docs.microsoft.com/en-us/cpp/standard-library/overloading-the-output-operator-for-your-own-classes?view=msvc-160
    https://www.csie.ntu.edu.tw/~sylee/courses/35500-99s/overload.htm

*/
#ifndef SIGNAL_HPP
#define SIGNAL_HPP

#include <complex>
#include <iostream>
#include <vector>

class Signal {
private:
    int origin_index;
    std::vector<int> vals;
    std::vector<std::vector<int>> get_matrix(const std::vector<int>& v1,
                                             const std::vector<int>& v2);
    std::vector<int> get_correlation_vals(
        const std::vector<std::vector<int>>& matrix);
    std::vector<std::vector<std::complex<double>>> multiply_matrices(
        const std::vector<std::vector<std::complex<double>>>& A,
        const std::vector<std::vector<int>>& B);

public:
    Signal();
    Signal(int origin_index, const std::vector<int>& vals);
    Signal shift(int k);
    Signal reverse();
    Signal reverse_and_shift(int k);
    // Signal time_scale(int a, int b);
    Signal time_scale(int c);
    Signal multiply_scalar(int scalar);
    // overloading *, +, and cout operators
    // int Signal::operator[] (int index);
    friend Signal operator*(const Signal& sig1, const Signal& sig2);
    friend Signal operator+(const Signal& sig1, const Signal& sig2);
    friend std::ostream& operator<<(std::ostream& os, const Signal& sig);

    Signal linear_convolution(const Signal& other);
    Signal circular_convolution(const Signal& other);
    Signal auto_correlate();
    Signal cross_correlate(const Signal& other);

    std::vector<std::complex<double>> DFT();
};

#endif