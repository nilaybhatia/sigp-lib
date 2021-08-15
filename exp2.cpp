#include <iostream>
#include <vector>
#include "Signal.hpp"

using namespace std;

int main(){
    cout << "Enter length of signal ";
    int n; cin >> n;
    cout << "Enter space-separated values of signal\n";
    vector<int> vals(n);
    for(int i = 0; i < n; i++) cin >> vals[i];
    cout << "Enter index of origin ";
    int origin; cin >> origin;
    Signal sig(origin, vals);

    cout << "\nShifted:\n";
    cout << "Enter value of k in x(+-n - k)";
    int k; cin >> k;
    Signal shifted = sig.shift(k);
    cout << shifted;
    cout << "Enter value of k ";
    cin >> k;
    shifted = sig.shift(k);
    cout << shifted;

    cout << "\nReversed:\n";
    cout << sig.reverse();

    cout << "\nReversed and shifted:\n";
    cout << "Enter value of k in x(+-n - k)";
    cin >> k;
    Signal reverse_shifted = sig.reverse_and_shift(k);
    cout << reverse_shifted;

    cout << "Enter value of k in x(+-n - k)";
    cin >> k;
    reverse_shifted = sig.reverse_and_shift(k);
    cout << reverse_shifted;

    // cout << "\nTime scaling:\n";
    // cout << "Enter the value of c in x(c * n)";
    // int c; cin >> c;
    // Signal time_scaled = sig.time_scale(c);
    // cout << time_scaled;

    // cout << "Enter the value of a, b in x(c * n)";
    // cin >> c;
    // time_scaled = sig.time_scale(c);
    // cout << time_scaled;


    cout << "\nScalar multiplication:\n";
    cout << "Enter scalar ";
    int scalar; cin >> scalar;
    Signal scaled = sig.multiply_scalar(scalar);
    cout << scaled;

    {
    cout << "\nSignal multiplication:\n";
    
    cout << "Enter length of signal 1";
    int l1; cin >> l1;
    cout << "Enter space-separated values of signal 1\n";
    vector<int> vals1(l1);
    for(int i = 0; i < l1; i++) cin >> vals1[i];
    cout << "Enter index of origin 1 ";
    int origin1; cin >> origin1;
    Signal sig1(origin1, vals1);
    
    cout << "Enter length of signal 2 ";
    int l2; cin >> l2;
    cout << "Enter space-separated values of signal 2\n";
    vector<int> vals2(l2);
    for(int i = 0; i < l2; i++) cin >> vals2[i];
    cout << "Enter index of origin 2";
    int origin2; cin >> origin2;
    Signal sig2(origin2, vals2);

    cout << sig1*sig2;
    }

    {
    cout << "\nSignal addition\n";
    cout << "Enter length of signal 1";
    int l1; cin >> l1;
    cout << "Enter space-separated values of signal 1\n";
    vector<int> vals1(l1);
    for(int i = 0; i < l1; i++) cin >> vals1[i];
    cout << "Enter index of origin 1 ";
    int origin1; cin >> origin1;
    Signal sig1(origin1, vals1);
    
    cout << "Enter length of signal 2 ";
    int l2; cin >> l2;
    cout << "Enter space-separated values of signal 2\n";
    vector<int> vals2(l2);
    for(int i = 0; i < l2; i++) cin >> vals2[i];
    cout << "Enter index of origin 2";
    int origin2; cin >> origin2;
    Signal sig2(origin2, vals2);
    
    cout << sig1 + sig2;
    
    }

        
    return 0;
}