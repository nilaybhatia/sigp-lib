#include "Signal.hpp"
using namespace std;

int main() {
    Signal<int> si(0, {1, 2, 3, 4});
    Signal<double> sd(0.0, {1.5, 2.5, 3.5, 4.5});
    return 0;
}