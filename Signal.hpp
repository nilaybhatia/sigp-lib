/*
Author: Nilay Bhatia (@nilaybhatia)
Date: 13th August 2021
References:
    https://docs.microsoft.com/en-us/cpp/cpp/header-files-cpp?view=msvc-160
    https://docs.microsoft.com/en-us/cpp/standard-library/overloading-the-output-operator-for-your-own-classes?view=msvc-160
    https://www.csie.ntu.edu.tw/~sylee/courses/35500-99s/overload.htm
    https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates

*/
#ifndef SIGNAL_HPP
#define SIGNAL_HPP

#include <complex>
#include <iostream>
#include <vector>

// Todo: must be templated
template <class T>
class Signal {
private:
    T origin_index;
    std::vector<T> vals;
    std::vector<std::vector<T>> get_matrix(const std::vector<T>& v1,
                                           const std::vector<T>& v2);
    std::vector<T> get_correlation_vals(
        const std::vector<std::vector<T>>& matrix);
    std::vector<std::vector<T>> multiply_matrices(
        const std::vector<std::vector<T>>& A,
        const std::vector<std::vector<T>>& B);

public:
    Signal();
    Signal(T origin_index, const std::vector<T>& vals);
    Signal<T> shift(int k);
    Signal<T> reverse();
    Signal<T> reverse_and_shift(int k);
    Signal<T> time_scale(int a, int b);
    Signal<T> time_scale(int c);
    Signal<T> multiply_scalar(int scalar);
    // overloading *, +, and cout operators
    // int Signal::operator[] (int index);
    friend Signal<T> operator*(const Signal<T>& sig1, const Signal<T>& sig2);
    friend Signal<T> operator+(const Signal<T>& sig1, const Signal<T>& sig2);
    friend std::ostream& operator<<(std::ostream& os, const Signal<T>& sig);

    Signal<T> linear_convolution(const Signal& other);
    Signal<T> circular_convolution(const Signal& other);
    Signal<T> auto_correlate();
    Signal<T> cross_correlate(const Signal& other);

    Signal DFT();
};

#endif