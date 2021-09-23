#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Signal.hpp"
using namespace std;

unsigned int lg(unsigned int num) {
    unsigned int ans = 0;
    while (num > 0) {
        ans++;
        num >>= 1;
    }
    return ans;
}
int main() {
    cout << fixed << setprecision(3);

    cout << "INPUT: \n";
    // 2-pt signal {1, 3}. Complex parts are 0
    ComplexSignal sg0(vector<complex<double>>({{1, 0}, {3, 0}}));
    cout << sg0;
    ComplexSignal fft0 = sg0.FFT();
    cout << "OUTPUT:\nIts FFT is: \n";
    int N0 = 2;
    std::cout << "Total real multiplications = " << 2 * N0 * lg(N0) << "\n";
    std::cout << "Total real additions = " << 3 * N0 * lg(N0) << "\n";
    cout << fft0 << "\n";

    cout << "INPUT: \n";
    ComplexSignal sg1(
        vector<complex<double>>({{1, 0}, {2, 0}, {3, 0}, {4, 0}}));
    cout << sg1;
    ComplexSignal fft1 = sg1.FFT();
    cout << "OUTPUT:\nIts FFT is: \n";
    int N1 = 4;
    std::cout << "Total real multiplications = " << 2 * N1 * lg(N1) << "\n";
    std::cout << "Total real additions = " << 3 * N1 * lg(N1) << "\n";
    cout << fft1 << "\n";

    cout << "INPUT: \n";
    ComplexSignal sg2(
        vector<complex<double>>({{1, 0}, {-1, 0}, {1, 0}, {-1, 0}}));
    cout << sg2;
    ComplexSignal fft2 = sg2.FFT();
    cout << "OUTPUT:\nIts FFT is: \n";
    int N2 = 4;
    std::cout << "Total real multiplications = " << 2 * N2 * lg(N2) << "\n";
    std::cout << "Total real additions = " << 3 * N2 * lg(N2) << "\n";
    cout << fft2;

    cout << "INPUT: \n";
    ComplexSignal sg3(vector<complex<double>>(
        {{1, 0}, {2, 0}, {3, 0}, {4, 0}, {4, 0}, {3, 0}, {2, 0}, {1, 0}}));
    cout << sg3;
    ComplexSignal fft3 = sg3.FFT();
    cout << "OUTPUT:\nIts FFT is: \n";
    int N3 = 8;
    std::cout << "Total real multiplications = " << 2 * N3 * lg(N3) << "\n";
    std::cout << "Total real additions = " << 3 * N3 * lg(N3) << "\n";
    cout << fft3;

    return 0;
}