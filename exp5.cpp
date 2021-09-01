#include <iostream>
#include <vector>
#include <complex>

#include "Signal.hpp"
using namespace std;

int main() {
    // vector<vector<int>> v1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    // vector<vector<int>> v2 = {{1, 2}, {4, 5}, {7, 8}};

    // Signal sg(0, {1, 2, 3});
    // auto result = sg.multiply_matrices(v1, v2);
    // for(int i = 0; i < result.size(); i++){
    //     for(int j = 0; j < result[0].size(); j++){
    //         cout << result[i][j] << " ";
    //     }
    //     cout << "\n";
    // }
    complex<int> c(2, 3);
    complex<int> d(2*c);
    d = 2 * c;
    cout << real(d) << "j" << imag(d) << "\n";
    cout << d;
    // Signal<int> sg(0, {1, 2, 3, 4});
    // cout << sg.DFT();
    return 0;
}