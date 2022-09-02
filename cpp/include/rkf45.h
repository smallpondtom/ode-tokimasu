#ifndef __RKF45_H__
#define __RKF45_H__

#include <Eigen/Dense>
#include <functional>
#include <math.h>

template<unsigned int N>
class RKF45
{
  public:
    using VectorNd = Eigen::Matrix<double, N, 1>;
    using ArrayNd = Eigen::Array<double, N, 1>;
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

    RKF45(std::function<VectorNd(double, VectorNd)> f, double stepsize, Eigen::Array2d timespan,
                                 double AbsTol = 1e-10, double RelTol = 1e-8) 
      : ode(f), h(stepsize), tspan(timespan), atol(AbsTol), rtol(RelTol) {};
    
    void call(double& t, VectorNd& x)
    {
      while (true)
      {
        // Make sure that end point is included
        float dt = this->tspan(1) - t;
        if (dt < this->h)
          this->h = dt;

        VectorNd xn = run(t, x, this->h);
        double En = nerror(xn);
        
        if (En <= 1)
        {
          if (this->ct < 2)
          {
            this->h *= h110(En);
            this->ct++;
          }
          else
            this->h *= h211(En, this->h);
          
          // Save latest error and stepsize values
          this->Enm1 = En;
          this->hnm1 = this->h;

          // Update states and time
          x = xn;
          t += h;
          
          break;
        }
        else
          this->h *= h110(En);
      }
    }

    VectorNd run(double t, VectorNd x, double h)
    {
      k1 = h * ode(t, x);
      k2 = h * ode(t + h*1.0/4, x + k1*1.0/4);
      k3 = h * ode(t + h*3.0/8, x + k1*3.0/32 + k2*9.0/32);
      k4 = h * ode(t + h*12.0/13, x + k1*1932.0/2197 - k2*7200.0/2197 + k3*7296.0/2197);
      k5 = h * ode(t + h, x + k1*439.0/216.0 - k2*8.0 + k3*3680.0/513.0 - k4*845.0/4104.0);
      k6 = h * ode(t + h/2.0, x - k1*8.0/27.0 + k2*2.0 - k3*3544.0/2565.0 + k4*1859.0/4104.0 - k5*11.0/40.0);

      VectorNd xn = x + (16.0/135.0*k1 + 6656.0/12825.0*k3 + 28561.0/56430.0*k4 - 9.0/50.0*k5 + 2.0/55.0*k6);  // order of 5 (better)
      this->err = 1.0/360*k1 - 128.0/4275*k3 - 2197.0/75240*k4 + 1.0/50*k5 + 2.0/55*k6;
      return xn;
    }
    
    /** Public variables */
    double h;  // stepsize

  protected:
    double h211(double En, double hn)
    {
      double eta = sf * pow(eps / En, a) * pow(Enm1 / eps, b) * pow(hn / hnm1, c);
      // Clamp eta 
      return (eta < 0.01) ? 0.01 : (eta > 2) ? 2.0 : eta;
    }

    double h110(double En)
    {
      double eta = sf * pow(eps / En, 1.0/(p+1));
      return eta;
    }


  private:
    double nerror(VectorNd x)
    {
      // Compute normalized error
      ArrayNd tmp = this->err.abs() / (atol + x.array().abs() * rtol);
      return std::max(1e-10, tmp.maxCoeff());
    }
    
    /** Private variables */
    std::function<VectorNd(double, VectorNd)> ode;
    VectorNd k1, k2, k3, k4, k5, k6;
    const double atol, rtol;
    const int p = 4;  // order of numerical integration;
    ArrayNd err;  
    Eigen::Array2d tspan;

    // Store error values and stepsize values
    unsigned int ct = 0;  // counter
    double Enm1 = 1.0;  // save previous normalized error
    double hnm1 = 1.0;  // save previous stepsize

    // Parameters for the stepsize controller
    const float eps = 0.8;
    const float sf = 0.9;  // safety factor
    const float a = 0.25/(p+1), b = -0.25/(p+1), c = -0.25;
};

#endif
