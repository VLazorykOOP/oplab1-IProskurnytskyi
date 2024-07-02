// g++ -std=c++11 -o algorithm main.cpp
// ./algorithm

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

struct Data {
    std::vector<double> X, T, U;
};

Data loadData(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Unable to open file");

    Data data;
    double x, t, u;
    while (file >> x >> t >> u) {
        data.X.push_back(x);
        data.T.push_back(t);
        data.U.push_back(u);
    }
    return data;
}

double interpolate(double x, const std::vector<double> &X, const std::vector<double> &Y) {
    auto it = std::lower_bound(X.begin(), X.end(), x);
    if (it == X.end()) return Y.back();
    if (it == X.begin()) return Y.front();

    size_t idx = it - X.begin();
    double x1 = X[idx - 1], x2 = X[idx];
    double y1 = Y[idx - 1], y2 = Y[idx];
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

double Grs(double x, double y);
double Gold(double x, double y);
double Glr(double x, double y);
double Srz(double x, double y, double z);
double T(double x);
double U(double x);

Data dataX_1_1, dataX1_00, dataX00_1;

double fun(double x, double y, double z) {
    double srz = Srz(x, y, z);
    double grs1 = Grs(x, y);
    double grs2 = Grs(y, z);
    double grs3 = Grs(z, x);
    return x * grs1 + y * grs2 + 0.33 * x * y * grs3;
}

double Grs(double x, double y) {
    if (x != y) {
        return 0.1389 * Srz(x + y, Gold(x, y), Glr(x, x * y)) +
               1.8389 * Srz(x - y, Gold(x, y / 5), Glr(5 * x, x * y)) +
               0.83 * Srz(x - 0.9, Glr(y, x / 5), Gold(5 * y, y));
    }
    return 0;
}

double Gold(double x, double y) {
    if (y > x) return y / x;
    if (x > y || x == y) return x / x;
    return 0;
}

double Glr(double x, double y) {
    if (std::abs(x) < 1) return y;
    if (std::abs(x) >= 1 || std::abs(y) >= 0.1) return std::sqrt(x * x + y * y - 4);
    return 0;
}

double Srz(double x, double y, double z) {
    double t_x = T(x);
    double u_z = U(z);
    double t_y = T(y);
    double u_y = U(y);
    return (t_x + u_z) * x * y - (t_y + u_y) * (x - y);
}

double T(double x) {
    if (x < -1) return interpolate(1 / x, dataX_1_1.X, dataX_1_1.T);
    if (x > 1) return interpolate(1 / x, dataX1_00.X, dataX1_00.T);
    return interpolate(x, dataX00_1.X, dataX00_1.T);
}

double U(double x) {
    if (x < -1) return interpolate(1 / x, dataX_1_1.X, dataX_1_1.U);
    if (x > 1) return interpolate(1 / x, dataX1_00.X, dataX1_00.U);
    return interpolate(x, dataX00_1.X, dataX00_1.U);
}

int main() {
    try {
        dataX_1_1 = loadData("dat_X_1_1.dat");
        dataX1_00 = loadData("dat_X1_00.dat");
        dataX00_1 = loadData("dat_X00_1.dat");

        double x, y, z;
        std::cout << "Enter x, y, z: ";
        std::cin >> x >> y >> z;

        double result = fun(x, y, z);
        std::cout << "fun(" << x << ", " << y << ", " << z << ") = " << result << std::endl;
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
    return 0;
}
