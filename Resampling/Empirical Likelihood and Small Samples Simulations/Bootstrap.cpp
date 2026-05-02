#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main() {
    vector<int> vec = {10, 20, 30};
    vec.push_back(40);

    for (int i = 0; i < vec.size(); ++i) {
        cout << vec[i] << " ";
    }
}

vector<int> bootstraptest(vector<int> x, int B) {
  int n = x.size();
  vector<int> results(B); // Initialize vector of length B
  
  for(int b = 0; b < B; ++b) {
    double sum = 0.0;
    for(int i = 0; i < n; ++i) {
      int idx = floor(R::runif(0, n));  // sample index with replacement
      sum += x[idx];
    }
    results[b] = sum / n;
  }
  
  return results;
}