#ifndef random_matrix
#define random_matrix

#include <random>
#include <algorithm>

template <typename T>
void fill_random_uniform_3d(T* arr, size_t dim1, size_t dim2, size_t dim3, double low, double high) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(low, high);

    size_t total = dim1 * dim2 * dim3;
    std::generate(arr, arr + total, [&]() { return dist(gen); });
}


template <typename T>
void fill_random_uniform_1d(T begin, T end, double low, double high) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(low, high);

    std::generate(begin, end, [&]() { return dist(gen); });
}

#endif 
