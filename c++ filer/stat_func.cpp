#include <iostream>
#include <cmath>

// Function to calculate factorial
unsigned long long factorial(int n) {
    unsigned long long result = 1;
    for (int i = 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}

// Function to calculate combinations (n choose k)
unsigned long long combinations(int n, int k) {
    return factorial(n) / (factorial(k) * factorial(n - k));
}

// Function to calculate hypergeometric probability
double hypergeometricProbability(int N, int K, int n, int k) {
    unsigned long long numerator = combinations(K, k) * combinations(N - K, n - k);
    unsigned long long denominator = combinations(N, n);
    return static_cast<double>(numerator) / denominator;
}

int main() {
    int N = 50; // Total number of items
    int K = 10; // Number of success items in the population
    int n = 5;  // Number of items drawn
    int k = 2;  // Number of success items drawn

    double probability = hypergeometricProbability(N, K, n, k);
    std::cout << "Hypergeometric probability: " << probability << std::endl;

    return 0;
}

