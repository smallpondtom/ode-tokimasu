#include <Eigen/Core>
#include <Eigen/Dense>
#include "exampleODE.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include "rkf45.h"  // CHANGE HERE!!!
#include <utility>


int main(int argc, char *argv[]) 
{
  // Setup for numerical integration 
  float ti = 0, tf = 30;
  Eigen::Array2d tspan(ti, tf);
  double h = 0.001;

  // Choose test case
  TEST2* test = new TEST2();  // CHANGE HERE!!!
  Eigen::Vector2d x_n = test->x0;  // numerical states
  Eigen::Vector2d x_e = test->x0;  // exact solution states
  
  // Initialize solver
  RKF45<2> rkf45(test->numerical, h, tspan);  // CHANGE HERE!!!
  // Setup csv file to store data
  std::ofstream results("numerical_and_exact_results.csv");

  std::cout << "[DEBUG] Start Numerical Integration\n";
  double time = ti;
  while (time <= tf)
  {
    // Exact solution
    x_e = test->exact(time);

    { // Data output
      std::printf("[INFO] %8.3f : (N) ", time);
      for (int i = 0; i < x_n.size(); i++)
      {
        std::printf("%12.8e ", x_n(i));
      }
      std::cout << " (E) ";
      for (int i = 0; i < x_e.size(); i++)
      {
        std::printf("%12.8e ", x_e(i));
      }
      std::cout << '\n';

      // Save data into a CSV file
      results << time;
      for (int i = 0; i < x_n.size(); i++)
      {
        results << "," << x_n(i);
      }
      for (int i = 0; i < x_e.size(); i++)
      {
        results << "," << x_e(i);
      }
      results << '\n';
    }
    
    // Numerical Solution
    rkf45.call(time, x_n);  // CHANGE HERE!!!
  } 

  std::cout << "[DEBUG] Done\n";
  return 0; 
}
