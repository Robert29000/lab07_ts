//
//  main.cpp
//  lab07
//
//  Created by Роберт Артур Меликян on 12/12/2020.
//  Copyright © 2020 Роберт Артур Меликян. All rights reserved.
//

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include "TextTable.h"

std::vector<double> f_values;
std::vector<double> x_values;
std::vector<double> noise_values;
std::vector<double> filter_values;

const double x_min = 0;
const double x_max = M_PI;
const int K = 100;
const double a = 0.25;
const int L = 10;
const double P = 0.95;
const double E = 0.01;
const int r = 3;
const auto M = (r - 1)/2;

double f_k(double x_k){
    return sin(x_k) + 0.5;
}

double x_k(double k){
    return x_min + k*(x_max - x_min)/K;
}

std::string get_string_vec(std::vector<double> vec){
    std::string res = "[ ";
    for (auto i = 0; i < vec.size(); i++){
        res += std::to_string(vec[i]) + " ";
    }
    res += "]";
    return res;
}

double filter_k(size_t k, const std::vector<double>& a){
    if (k < M || k > K - M){
        return 0.0;
    }
    double res = 0.;
    for (auto j = k - M; j <= k + M; j++){
        res += a[j + M - k] / noise_values[j];
    }
    return pow(res, -1);
}

double w(const std::vector<double>& coefs){
    double res = 0;
    for (auto k = 1; k <= K; k++){
        res += fabs(filter_k(k, coefs) - filter_k(k - 1, coefs));
    }
    return res;
}

double d(const std::vector<double>& coefs){
    double res = 0;
    for (auto k = 0; k <= K; k++){
        res += fabs(filter_k(k, coefs) - noise_values[k]);
    }
    return 1.0/K * res;
}

double J(double y, const std::vector<double>& coefs){
    return y * w(coefs) + (1 - y) * d(coefs);
}

double sum_a(size_t from, size_t to, std::vector<double> a){
    double res = 0;
    for (auto i = from; i <= to; i++){
        res += a[i];
    }
    return res;
}


std::vector<double> get_a_coefs(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1);
    std::vector<double> a(r);
    a[M] = dis(gen);
    for (auto m = 1; m <= M; m++){
        std::uniform_real_distribution<double> dis(0, 1 - sum_a(m, r - m - 1, a));
        a[m - 1] = a[r - m] = 0.5 * dis(gen);
    }
    a[0]=a[r-1]=0.5 * (1 - sum_a(1, r-2, a));
    return a;
}

int main(int argc, const char * argv[]) {
    
    TextTable exp('-', '|', '+');
    TextTable res('-', '|', '+');
    exp.add("y");
    exp.add("dist");
    exp.add("alpha");
    exp.add("w");
    exp.add("d");
    exp.endOfRow();
    
    res.add("y*");
    res.add("J");
    res.add("w");
    res.add("d");
    res.endOfRow();
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-a, a);
    
    f_values.reserve(K+1);
    x_values.reserve(K+1);
    noise_values.reserve(K+1);
    filter_values.reserve(K+1);
    for (auto k = 0; k <= K; k++){
        auto x = x_k(k);
        auto f = f_k(x);
        auto n = f + dis(gen);
        x_values.push_back(x);
        f_values.push_back(f);
        noise_values.push_back(n);
    }
    auto N = ceil(log(1-P)/log(1 - (E/(x_max - x_min))));
    double total_dist_min = 10000;
    std::vector<double> total_coefs;
    double total_y = 0;
    double j_total = 0;
    for (auto l = 0; l <= L; l++){
        auto y = double(l)/L;
        std::vector<double> local_coefs;
        double j_local = 10000;
        for (auto i = 0; i < N; i++){
            auto a = get_a_coefs();
            auto j_val = J(y, a);
            if (j_val < j_local){
                j_local = j_val;
                local_coefs = a;
            }
        }
        double local_dist_min = fabs(w(local_coefs)) + fabs(d(local_coefs));
        exp.add(std::to_string(y));
        exp.add(std::to_string(local_dist_min));
        exp.add(get_string_vec(local_coefs));
        exp.add(std::to_string(w(local_coefs)));
        exp.add(std::to_string(d(local_coefs)));
        exp.endOfRow();
        if (local_dist_min < total_dist_min){
            total_dist_min = local_dist_min;
            total_coefs = local_coefs;
            total_y = y;
            j_total = j_local;
        }
    }
    
    res.add(std::to_string(total_y));
    res.add(std::to_string(j_total));
    res.add(std::to_string(w(total_coefs)));
    res.add(std::to_string(d(total_coefs)));
    res.endOfRow();
    
    std::cout << exp << std::endl;
    std::cout << res << std::endl;
    for (auto k = 0; k <= K; k++){
        std::cout << "(" << x_values[k] << ";" << noise_values[k] << ") ";
    }
    std::cout << std::endl;
    for (auto k = 0; k <= K; k++){
        std::cout << "(" << x_values[k] << ";" << filter_k(k, total_coefs) << ") ";
    }
    std::cout << std::endl;
    return 0;
}
