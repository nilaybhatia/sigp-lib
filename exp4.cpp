#include <iostream>
#include <vector>
#include "Signal.hpp"

using namespace std;

int main(){
    cout << "Auto correlation\n";
    Signal x(0, {1, 2, 1, 2});
    cout << x.auto_correlate();
    Signal s(1, {1, 2, 1, 2});
    cout << s.auto_correlate();
    // Conclusion: Shifting a signal does not affect
    // the auto-correlation

    cout << "Cross correlation\n";
    cout << "Larger signal must be followed by smaller signal";
    Signal x1(0, {1, 2, 3, 4});
    Signal x2(0, {5, 0, 6});
    cout << x1.cross_correlate(x2);

    Signal s1(0, {1, 2, 1, 2});
    Signal s2(0, {0, 1, 2, 1, 2}); // above signal right-shifted by 1
    cout << s2.cross_correlate(s1);
    // Conclusion: When shifting a signal towards right by k units,
    // and then taking cross-correaltion, the origin shifts by k units towards left
    // from the auto-correlated signal
    return 0;
}