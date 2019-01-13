//
//  Utils.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#include "Random.hpp"

using namespace std;

double randomQMC(int axis, long long i) {
    int base = prime[axis];
    double f = 1, res = 0;
    while (i > 0) {
        f /= base;
        res += f * (i % base);
        i /= base;
    }
    return res;
}

static double random01() {
    static mt19937 *generator = nullptr;
    if (!generator)
        generator = new mt19937((unsigned int)clock());
    static uniform_real_distribution<> dis(0, 1);
    return dis(*generator);
}

double randomdist(double l, double r, int axis, long long i) {
    if (axis == -1)
        return l + random01() * (r - l);
    return l + randomQMC(axis, i) * (r - l);
}
