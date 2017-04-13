#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if (estimations.size() == 0) {
    return rmse;
  }
  if (ground_truth.size() != estimations.size()) {
    return rmse;
  }

  vector<VectorXd> diff(estimations.size());
  
  //Calculate the difference between estimation and ground truth
  std::transform(estimations.begin(), estimations.end(), ground_truth.begin(), \
                diff.begin(), [](VectorXd d1, VectorXd d2) { return d1 - d2;});
                
  //Calculate the square of every element
  std::transform(diff.begin(), diff.end(), diff.begin(), 
                [](VectorXd d){return (VectorXd)(d.array()*d.array());});
  
  // Calculate the sum of all the element
  rmse = std::accumulate(diff.begin(), diff.end(), rmse, 
                [](VectorXd d1, VectorXd d2) { return d1 + d2;} );

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj = MatrixXd::Zero(3,4);
  
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];
  
  float pxy2= px*px + py*py;
  
  float pxy2_root = sqrt(pxy2);
  float pxy2_3_root = pow(pxy2_root, 3);
  
  Hj(0,0) = px/pxy2_root;
  Hj(0,1) = py/pxy2_root;
  Hj(1,0) = - py/pxy2;
  Hj(1,1) = px/pxy2;
  Hj(2,0) = py*(vx * py - vy * px)/pxy2_3_root;
  Hj(2,1) = px*(vy * px - vx * py)/pxy2_3_root;
  Hj(2,2) = px/pxy2_root;
  Hj(2,3) = py/pxy2_root;
  
  return Hj;
}
