#ifndef __RK4_H__
#define __RK4_H__

#include <Eigen/Dense>
#include <functional>

template<unsigned int N>
class RK4
{
  public:
    using VectorNd = Eigen::Matrix<double, N, 1>;

    RK4(std::function<VectorNd(double, VectorNd)> f) : ode(f) {}

    void call(double& t, VectorNd& x, double h)
    {
      k1 = h * ode(t, x);
      k2 = h * ode(t + h / 2, x + k1 / 2);
      k3 = h * ode(t + h / 2, x + k2 / 2);
      k4 = h * ode(t + h, x + k3);

      x += (k1 + k2 * 2 + k3 * 2 + k4) / 6;
      t += h;
      return;
    }

  private:
    std::function<VectorNd(double, VectorNd)> ode;
    VectorNd k1, k2, k3, k4;
};

// THE LINEAR ALGEBRA METHOD DID NOT WORK OUT FOR. NEED MORE PRACTICE WITH EIGEN....

/* HEADER FILE
Eigen::VectorXd rk4(Eigen::VectorXd (*derivative)(double, Eigen::VectorXd),
                    double t, double h, Eigen::VectorXd x);

const double Aval[] = {0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0};
const Eigen::Array44d A(Aval);

const double Bval[] = {1.0/6, 1.0/3, 1.0/3, 1.0/6};
const Eigen::Array4d B(Bval);

const double Cval[] = {0, 0.5, 0.5, 1};
const Eigen::Array4d C(Cval);
*/

/* SOURCE FILE
#include <Eigen/Core>
#include <Eigen/Dense>
#include "rk4.h"
#include <iostream>

Eigen::VectorXd rk4(Eigen::VectorXd (*derivative)(double, Eigen::VectorXd),
                    double t, double h, Eigen::VectorXd x)
{
  int order = 4;
  int N = x.size();
  Eigen::ArrayXXd K = Eigen::MatrixXd::Zero(order, N);

  for (int i = 0; i < order; i++)
  {
    Eigen::VectorXd Arow = A.row(i);
    Eigen::VectorXd tmp = derivative(t + h*C(i),
        x.array() + h * (K.array().colwise() * Arow.array()).colwise().sum()
    );
  }
 
  Eigen::VectorXd xnew;
  xnew = x.array() + h * (K.array().colwise() * B.array()).colwise().sum();

  return xnew;
}
*/

#endif
