#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include "Signal.hpp"

Signal::Signal(){
    *this = Signal(0, std::vector<int>{});
}

Signal::Signal(int origin_index, const std::vector<int>& vals) {
    this->origin_index = origin_index;
    this->vals = vals;
}

Signal Signal::shift(int k){ 
    // std::cout << *this << "k = " << k;
    std::deque<int> new_seq((this->vals).begin(), (this->vals).end());
    int new_origin_index;
    if(k > 0){
        new_origin_index = std::max(this->origin_index - k, 0);
        for(int i = 1; i <= k - this->origin_index; i++) new_seq.push_front(0);
    }
    else{
        new_origin_index = this->origin_index - k;
        for(int i = 1; i <= new_origin_index - (int)this->vals.size() + 1; i++) new_seq.push_back(0);
    }
    return Signal(new_origin_index, std::vector<int>(new_seq.begin(), new_seq.end()));
}

Signal Signal::reverse(){
    std::vector<int> reversed(this->vals.rbegin(), this->vals.rend());
    int new_origin_index = this->vals.size() - 1 - this->origin_index;
    return Signal(new_origin_index, reversed);
}

Signal Signal::reverse_and_shift(int k){
    Signal reversed = this->reverse();
    Signal reverse_shifted = reversed.shift(-k);
    return reverse_shifted;
}

Signal Signal::time_scale(int c){
    // TODO: Support any rational like a/b
    int new_origin_index = this->origin_index;
    std::deque<int> new_vals;
    new_vals.push_back(this->vals[new_origin_index]);
    
    // Shrinks
    for(int i = new_origin_index + 1; i < this->vals.size(); i++){
        if((i - new_origin_index) % c == 0) new_vals.push_back(this->vals[i]);
    }
    for(int i = new_origin_index - 1; i >= 0; i++){
        if((new_origin_index - i) % c == 0) new_vals.push_front(this->vals[i]);
    }
    return Signal(new_origin_index, std::vector<int> (new_vals.begin(), new_vals.end()));
}

Signal Signal::multiply_scalar(int scalar){
    std::vector<int> new_vals((this->vals).begin(), (this->vals).end());
    for(int i = 0; i < new_vals.size(); i++){
        new_vals[i] *= scalar;
    }
    return Signal(this->origin_index, new_vals);
}

// int Signal::operator[] (int index){
//     if(index < 0 or index >= this->vals.size()) return 0;
//     else return this->vals[index];
// }

Signal operator* (const Signal& sig1, const Signal& sig2){
    int new_orgin_index = std::min(sig1.origin_index, sig2.origin_index);
    std::deque<int> new_vals;
    int i = sig1.origin_index, j = sig2.origin_index;

    while(i >= 0 && j >= 0){
        new_vals.push_front(sig1.vals[i] * sig2.vals[j]);
        i--;
        j--;
    }

    i = sig1.origin_index + 1, j = sig2.origin_index + 1;
    while(i < sig1.vals.size() && j < sig2.vals.size()){
        new_vals.push_back(sig1.vals[i] * sig2.vals[j]);
        i++;
        j++;
    }
    return Signal(new_orgin_index, std::vector<int>{new_vals.begin(), new_vals.end()});
}

Signal operator+ (const Signal& sig1, const Signal& sig2){
    int new_orgin_index = std::max(sig1.origin_index, sig2.origin_index);
    std::deque<int> new_vals;
    int i = sig1.origin_index, j = sig2.origin_index;

    while(i >= 0 || j >= 0){
        int a = (i < 0? 0 : sig1.vals[i]);
        int b = (j < 0? 0 : sig2.vals[j]);
        new_vals.push_front(a + b);
        i--;
        j--;
    }

    i = sig1.origin_index + 1, j = sig2.origin_index + 1;
    while(i < sig1.vals.size() || j < sig2.vals.size()){
        int a = (i >= sig1.vals.size()? 0 : sig1.vals[i]);
        int b = (j >= sig2.vals.size()? 0 : sig2.vals[j]);
        new_vals.push_back(a + b);
        i++;
        j++;
    }
    return Signal(new_orgin_index, std::vector<int>(new_vals.begin(), new_vals.end()));
}


std::ostream& operator<< (std::ostream& os, const Signal& seq){
    os << "The signal is \n";
    for(int ele : seq.vals) os << ele << " "; os << "\t";
    os << "with origin as " << seq.origin_index << "\n";
    return os;
}


Signal Signal::linear_convolution(const Signal& other){
    int new_origin_index = this->origin_index + other.origin_index;
    int n = this->vals.size(), m = other.vals.size();

    // n x m matrix
    std::vector<std::vector<int>> matrix(n, std::vector<int>(m));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            matrix[i][j] = this->vals[i] * other.vals[j];
        }
    }

    std::vector<int> result;
    for(int i = 0; i < n; i++){
        int sum = 0;
        for(int row = i, col = 0; row >= 0; row--, col++){
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    for(int j = 1; j < m; j++){
        int sum = 0;
        for(int row = n-1, col = j; col < m; row--, col++){
            sum += matrix[row][col];
        }
        result.push_back(sum);
    }
    Signal convoluted(new_origin_index, result);
    return convoluted;
}

