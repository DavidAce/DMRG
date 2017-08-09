//
// Created by david on 4/19/17.
//

#ifndef TRAINING_FUNCS_H
#define TRAINING_FUNCS_H


#include <iostream>
#include <iterator>
#include <vector>

template<typename T1, typename T2>
inline  long    mod (const T1 x, const T2 y){
    return x >= 0 ? x%y : x%y + y;
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    if ( !v.empty() ) {
        out << "[ ";
        std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

/** x^p for integers x and p using recursion: */
template<typename IntegerType>
IntegerType ipow(IntegerType x, IntegerType p){
    if (p == 0) return 1;
    if (p == 1) return x;
    int tmp = ipow(x, p/2);
    if (p%2 == 0) return tmp * tmp;
    else return x * tmp * tmp;
}





//
//template<typename Input, typename From, typename To>
//auto prod(const Input &in, const From from, const To to){
//    return accumulate(in.data() + from, in.data()+to,1,multiplies<>());
//};


#endif //TRAINING_FUNCS_H
