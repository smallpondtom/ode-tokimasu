#ifndef __DOP54_H__
#define __DOP54_H__

#include <Eigen/Dense>
#include <functional>
#include <math.h>
#include <algorithm>

template<unsigned int N>
class DOP54
{
  public:
    using VectorNd = Eigen::Matrix<double, N, 1>;
    using ArrayNd  = Eigen::Array<double, N, 1>;

    // Constructor: takes the ODE function, an initial stepsize, a timespan [t0, tf],
    // and optional absolute and relative tolerances.
    DOP54(std::function<VectorNd(double, VectorNd)> f, double stepsize,
          Eigen::Array2d timespan, double AbsTol = 1e-10, double RelTol = 1e-8)
      : ode(f), h(stepsize), tspan(timespan), atol(AbsTol), rtol(RelTol) {}

    // Adaptive integration step:
    // Advances the solution from time t by h (adjusting h as needed) while computing the error.
    // On a successful step, both the state x and the current time t are updated.
    void call(double &t, VectorNd &x)
    {
      while (true)
      {
        // Ensure we don't step past the final time.
        double dt = tspan(1) - t;
        if (dt < h)
          h = dt;

        // Compute one integration step and the embedded error estimate.
        VectorNd x_new = run(t, x, h);
        double En = nerror(x_new);

        // If the normalized error is acceptable, accept the step and update.
        if (En <= 1.0)
        {
          if (ct < 2)
          {
            h *= h110(En);
            ct++;
          }
          else
          {
            h *= h211(En, h);
          }
          // Save latest error and step size values (for use in h211)
          Enm1 = En;
          hnm1 = h;
          x = x_new;
          t += h;
          break;
        }
        else
        {
          // Reject step and reduce h.
          h *= h110(En);
        }
      }
    }

    // Run one step of size h: computes the 7 stages, updates the state with the 5th order estimate,
    // and computes an error estimate.
    VectorNd run(double t, VectorNd x, double h)
    {
      k1 = h * ode(t, x);
      k2 = h * ode(t + h * 1.0/5, x + k1 * 1.0/5);
      k3 = h * ode(t + h * 3.0/10, x + k1 * 3.0/40 + k2 * 9.0/40);
      k4 = h * ode(t + h * 4.0/5, x + k1 * 44.0/45 - k2 * 56.0/15 + k3 * 32.0/9);
      k5 = h * ode(t + h * 8.0/9, x + k1 * 19372.0/6561 - k2 * 25360.0/2187 +
                                  k3 * 64448.0/6561 - k4 * 212.0/729);
      k6 = h * ode(t + h, x + k1 * 9017.0/3168 - k2 * 355.0/33 +
                                  k3 * 46732.0/5247 + k4 * 49.0/176 - k5 * 5103.0/18656);
      k7 = h * ode(t + h, x + k1 * 35.0/384 + k3 * 500.0/1113 +
                                  k4 * 125.0/192 - k5 * 2187.0/6784 + k6 * 11.0/84);

      // 5th order solution estimate.
      VectorNd x_new = x + k1 * 5179.0/57600 + k3 * 7571.0/16695 + k4 * 393.0/640 -
                        k5 * 92097.0/339200 + k6 * 187.0/2100 + k7 * 1.0/40;
      // Embedded error estimate.
      err = k1 * (-71.0/57600) + k3 * (71.0/16695) - k4 * (71.0/1920) +
            k5 * (17253.0/339200) - k6 * (22.0/525) + k7 * (1.0/40);
      return x_new;
    }

    // (Optional) Retrieve the most recent error estimate vector.
    VectorNd getError() const { return err; }

  protected:
    // Step size adjustment factors.
    double h211(double En, double hn)
    {
      double eta = sf * pow(eps / En, a) * pow(Enm1 / eps, b) * pow(hn / hnm1, c);
      return (eta < 0.01) ? 0.01 : (eta > 2.0) ? 2.0 : eta;
    }

    double h110(double En)
    {
      double eta = sf * pow(eps / En, 1.0/(p+1));
      return eta;
    }

  private:
    // Compute the normalized error based on the embedded error estimate.
    double nerror(const VectorNd& x_new)
    {
      ArrayNd tmp = err.array().abs() / (atol + x_new.array().abs() * rtol);
      return std::max(1e-10, tmp.maxCoeff());
    }

    std::function<VectorNd(double, VectorNd)> ode;
    VectorNd k1, k2, k3, k4, k5, k6, k7;
    VectorNd err;
    double h;                     // Current stepsize.
    Eigen::Array2d tspan;         // Integration interval [t0, tf].
    const double atol, rtol;      // Error tolerances.
    const int p = 4;              // Order used for error control.

    // Variables for step size control history.
    unsigned int ct = 0;
    double Enm1 = 1.0;
    double hnm1 = 1.0;
    const float eps = 0.8;
    const float sf = 0.9;
    const float a = 0.25/(p+1);
    const float b = -0.25/(p+1);
    const float c = -0.25;
};

#endif
