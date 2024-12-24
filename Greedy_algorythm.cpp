#include <iostream>
#include <vector>
#include <algorithm>


int maxSatisfaction(std::vector<int>& satisfaction) {
    // методы rbegin и rend позволяют отсортировать в обратном порядке
    std::sort(satisfaction.rbegin(), satisfaction.rend());

    int total = 0;
    long result = 0;
    
    int n = size(satisfaction);
    // На каждом шаге проверяем, не становится ли total < 0
    for (size_t i = 0; i < n; i++){

        if (total + satisfaction[i] < 0){
            break;
        }
        
        total += satisfaction[i];
        result += total;
    }

    return result;
}


int main(){

    std::vector<int> inputArray = {-8,-1,0,5,-9};

    std::cout << maxSatisfaction(inputArray) << '\n';

    return 0;
}