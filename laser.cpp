#include "params.h"
#include "solver.h"

#include <cmath>

#if defined(PLANE)
double Simu::get_inject(double x, double y, double t)
{
    double omega = 0.2 * M_PI;
    return std::sin(omega * (x - t));
}
#elif defined(GAUSSIAN)
double Simu::get_inject(double x, double y, double t)
{
    const double a0 = 30.0;        // 场强幅度
    const double w0 = 40.0;        // 束腰半径
    const double focus = 200.0;    // 焦点位置
    const double l0 = 5.0;         // 波长
    const double FWHM = 100.0;     // 脉冲持续时间
    const double y_center = 150.0; // 束流横向中心位置
    const double t0 = 200.0;
    const double k0 = 2 * M_PI / l0;
    const double omega0 = k0; // 角频率（光速c=1）

    // 光束参数计算
    const double x_R = M_PI * w0 * w0 / l0;              // 瑞利长度
    const double r2 = pow(y - y_center, 2);              // 横向距离平方
    const double x_ = x - focus;                         // 离焦量
    const double wx = w0 * sqrt(1.0 + pow(x_ / x_R, 2)); // 光束半径
    const double Rx = x_ * (1.0 + pow(x_R / x_, 2));     // 波前曲率半径
    const double gouy = -atan(x_ / x_R);                 // Gouy相位

    // 空间包络
    const double spatial = exp(-r2 / (wx * wx));

    // 相位项
    const double phase = k0 * x - omega0 * t  // 基础相位
                         + k0 * r2 / (2 * Rx) // 波前曲率相位
                         + gouy;              // Gouy相位

    // 时间包络（含传播延迟修正）
    // const double delay = z - k0 * r2 / (2 * omega0 * Rx) + t0; // 总延迟（传播+波前曲率）
    const double delay = t0 + x_ + r2 / (2 * Rx); // 总延迟（传播+波前曲率）
    const double temporal = exp(-4 * log(2) * pow((t - delay) / FWHM, 2));

    return a0 * w0 / wx * spatial * temporal * cos(phase);
    // return a0 * w0 / wz * spatial * cos(phase);
}
#elif defined(GAUSSIAN_2)
double Simu::get_inject(double x, double y, double t)
{
    const double w0 = 40.0;
    const double tau = 20.0;
    const double a0 = 30.0;
    const double focus = 200;
    const double l0 = 5.0;
    const double k0 = 2 * M_PI / l0;

    const double xR = M_PI * w0 * w0 / l0;
    const double r2 = (y - 150) * (y - 150);
    const double x_ = focus - x;
    const double wx = w0 * std::sqrt(1.0 + (x_ / xR) * (x_ / xR));
    const double Rx = x_ * (1.0 + (xR / x_) * (xR / x_));
    const double gouy = -std::atan(x_ / xR);
    const double phi = k0 * r2 / (2.0 * Rx) + gouy + t ;
    const double prof =
        std::exp(4.0 * std::log(2.0) *
                 (-std::pow(((t + k0 * r2 / (2.0 * Rx)) / (2.0 * M_PI) - std::sqrt(2.0) * tau), 2) / (tau * tau)));

    return a0 * w0 / wx * std::exp(-r2 / (wx * wx)) * std::cos(phi) * prof;
}
#endif
