#include <iostream>
#include <vector>
#include "Signal.hpp"

using namespace std;
int main(){
    cout << "Enter number of elements in x(n)";
    int len_x; cin >> len_x;
    cout << "Enter x(n)";
    vector<int> vals_x(len_x);
    for(int i = 0; i < len_x; i++) cin >> vals_x[i];
    cout << "Enter origin"; 
    int origin_x; cin >> origin_x;
    Signal x(origin_x, vals_x);

    cout << "Enter number of elements in h(n)";
    int len_h; cin >> len_h;
    cout << "Enter h(n)";
    vector<int> vals_h(len_h);
    for(int i = 0; i < len_h; i++) cin >> vals_h[i];
    cout << "Enter origin"; 
    int origin_h; cin >> origin_h;
    Signal h(origin_h, vals_h);

    Signal result = x.linear_convolution(h);
    cout << result;
}