//#pragma once

#include <iostream>
#include <vector>

double f(double t, double y, double z, double a, double b, double c, double d) {
    return (d - b*z - c*y) / a;
}

std::vector <double> runge_kutta_4th_order_servo_motor(double ei, double La, double Ra, double Jl, double Jm, double Kt, double Ke, double y0, double dy0, double t0, double T){
    double a = (La / Ra) * (Jm / Kt);
    double b = Jm / Kt;
    double c = Ke;
    double d = ei;

    double Time = T;
    double tau = 0.1;
    double n = Time / tau;

    std::vector<double> y(n + 1), z(n + 1), t(n + 1);
    y[0] = y0;
    z[0] = dy0;
    t[0] = t0;

    for(size_t i = 0; i < n; i++){
        double k1_y = tau * z[i];
        double k1_z = tau * f(t0 + i * tau, y[i], z[i], a, b, c, d);

        double k2_y = tau * (z[i] + k1_z / 2);
        double k2_z = tau * f(t0 + i * tau + tau/2, y[i] + k1_y / 2, z[i] + k1_z / 2, a, b, c, d);

        double k3_y = tau * (z[i] + k2_z / 2);
        double k3_z = tau * f(t0 + i * tau + tau / 2, y[i] + k2_y / 2, z[i] + k2_z / 2, a, b, c, d);

        double k4_y = tau * (z[i] + k3_z);
        double k4_z = tau * f(t0 + i * tau + tau, y[i] + k3_y, z[i] + k3_z, a, b, c, d);

        y[i + 1] = y[i] + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
        z[i + 1] = z[i] + (k1_z + 2 * k2_z + 2 * k3_z + k4_z) / 6;
        t[i + 1] = t[i] + tau;
    }

    for(size_t i = 0; i < n; i++){
        std::cout << t[i] << " " << y[i] << std::endl;
    }
    return y;
}

int main(){
    runge_kutta_4th_order_servo_motor(1, 1, 1, 1, 1, 1, 1, 1, 10, 0, 5);
}




