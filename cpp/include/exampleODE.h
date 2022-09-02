#ifndef EXAMPLEODE_HPP
#define EXAMPLEODE_HPP

#include <Eigen/Dense>
#include <math.h>

class TEST1
{
  public:
    using Vector2d = Eigen::Matrix<double, 2, 1>;

    Vector2d x0 = { 1, 0.5 };  // initial conditions

    static Vector2d numerical(double t, Vector2d y)
    {
      Vector2d ynew(2);
      ynew(0) = y(1);
      ynew(1) = -4.0 * t * y(1) - (2 + 4*pow(t, 2))*y(0);
      return ynew;
    }

    static Vector2d exact(double t)
    {
      Vector2d ynew(2);
      ynew(0) = 0.5*exp(-pow(t, 2))*(t+2);
      ynew(1) = (-pow(t, 2) - 2*t + 0.5)*exp(-pow(t, 2));
      return ynew;
    }
};


class TEST2
{
  public:
    using Vector2d = Eigen::Matrix<double, 2, 1>;

    Vector2d x0 = { 2, 1 };  // initial conditions

    static Vector2d numerical(double t, Vector2d y)
    {
      Vector2d ynew(2);
      ynew(0) = (cos(t) - sin(t) * y(0)) / y(1);
      ynew(1) = sin(t);
      return ynew;
    }

    static Vector2d exact(double t)
    {
      Vector2d ynew(2);
      ynew(0) = (sin(t) + 2.0) / (-cos(t) + 2.0);
      ynew(1) = -cos(t) + 2.0;
      return ynew;
    }
};

#endif
