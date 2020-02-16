#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd ans(4);
  ans << 0,0,0,0;
  
  if(estimations.size() == 0){
    std::cout << "ERROR - The estimations vector is empty" << std::endl;
    return ans;
  }

  if(ground_truth.size() == 0){
    std::cout << "ERROR - The ground-truth vector is empty" << std::endl;
    return ans;
  }

  unsigned int n = estimations.size();
  
  if(n != ground_truth.size()){
    std::cout << "ERROR - The ground-truth and estimations vectors must have the same size." << std::endl;
    return ans;
  }
  
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array()*diff.array();
    ans += diff;
  }

  ans = ans / n;
  ans = ans.array().sqrt();
  return ans;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double c1 = px*px+py*py;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);
  
  if(fabs(c1) < 0.0001){
      std::cout << "ERROR - Division by Zero" << std::endl;
      return Hj;
  }
  
  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
