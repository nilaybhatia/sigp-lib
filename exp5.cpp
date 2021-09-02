#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Signal.hpp"
using namespace std;

int main() {
    cout << fixed << setprecision(3);
    Signal sg1(0, {1, 2, 3, 4});
    cout << sg1;
    vector<complex<double>> dft1 = sg1.DFT();
    cout << "It's DFT is: \n";
    for (auto ele : dft1)
        cout << noshowpos << real(ele) << showpos << " " << imag(ele) << "j"
             << noshowpos << "\t"
             << "(Amplitude: " << abs(ele) << ", Phase: " << arg(ele) << ")\n";
    cout << "\n";

    Signal sg2(0, {1, 2, 3, 4, 0, 0, 0, 0});
    cout << sg2;
    vector<complex<double>> dft2 = sg2.DFT();
    cout << "It's DFT is: \n";

    for (auto ele : dft2)
        cout << noshowpos << real(ele) << showpos << " " << imag(ele) << "j"
             << noshowpos << "\t"
             << "(Amplitude: " << abs(ele) << ", Phase: " << arg(ele) << ")\n";
    cout << "\n";
    // Alternate values are same, more values get inserted for greater
    // precision.

    Signal sg3(0, {1, 2, 3, 4, 4, 3, 2, 1});
    cout << sg3;
    vector<complex<double>> dft3 = sg3.DFT();
    cout << "It's DFT is: \n";

    for (auto ele : dft3)
        cout << noshowpos << real(ele) << showpos << " " << imag(ele) << "j"
             << noshowpos << "\t"
             << "(Amplitude: " << abs(ele) << ", Phase: " << arg(ele) << ")\n";
    cout << "\n";

    return 0;
}