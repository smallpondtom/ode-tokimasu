#ifndef __DOP54_H__
#define __DOP54_H__

#include <Eigen/Dense>
#include <functional>

template<unsigned int N>
class DOP54
{
  public:
    using VectorNd = Eigen::Matrix<double, N, 1>;
    using Vector7d = Eigen::Matrix<double, 7, 1>;
    using Matrix7d = Eigen::Matrix<double, 7, 7>;

    DOP54(std::function<VectorNd(double, VectorNd)> f) : ode(f) {};

    VectorNd call(double t, VectorNd& x, double h)
    {
      k1 = h * ode(t, x);
      k2 = h * ode(t + h*1.0/5, x + k1*1.0/5);
      k3 = h * ode(t + h*3.0/10, x + k1*3.0/40 + k2*9.0/40);
      k4 = h * ode(t + h*4.0/5, x + k1*44.0/45 - k2*56.0/15 + k3*32.0/9);
      k5 = h * ode(t + h*8.0/9, x + k1*19372.0/6561 - k2*25360.0/2187 + k3*64448.0/6561 - k4*212.0/729);
      k6 = h * ode(t + h, x + k1*9017.0/3168 - k2*355.0/33 + k3*46732.0/5247 + k4*49.0/176 - k5*5103.0/18656);
      k7 = h * ode(t + h, x + k1*35.0/384 + k3*500.0/1113 + k4*125.0/192 - k5*2187.0/6784 + k6*11.0/84);

      x += k1*5179.0/57600 + k3*7571.0/16695 + k4*393/640 - k5*92097.0/339200 + k6*187.0/2100 + k7*1/40;
      VectorNd err = k1 * -71.0/57600 + k3*71.0/16695 - k4*71.0/1920 + k5*17253.0/339200 - k6*22.0/525 + k7*1.0/40;
      return err;
    }

  private:
    std::function<VectorNd(double, VectorNd)> ode;
    VectorNd k1, k2, k3, k4, k5, k6, k7;
};

#endif

