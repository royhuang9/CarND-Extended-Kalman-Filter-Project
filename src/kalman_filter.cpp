#include "kalman_filter.h"
#include "tools.h"
#include <cmath>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose()*S.inverse();

  
  x_ = x_ + K*y;
  
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  Tools tools;
  MatrixXd Hj = tools.CalculateJacobian(x_);
  
  VectorXd y;
  VectorXd h_prime(3);
  
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  
  double pxy2 = px*px + py*py;
  
  double pxy2_root = sqrt(pxy2);
  
  h_prime[0] = pxy2_root;
  h_prime[1] = atan2(py, px);
  h_prime[2] = (px * vx + py * vy) / pxy2_root;
  y = z - h_prime;
  
  while (y(1) < -M_PI) {
    y(1) += 2*M_PI;
  }

  while (y(1) > M_PI) {
    y(1) -= 2*M_PI;
  }

  //cout<<"y(1):"<<y(1)<<",h_prime(1):"<<h_prime(1)<<"\n"<<endl;
  
  MatrixXd S = Hj * P_ * Hj.transpose() + R_;
  MatrixXd K = P_ * Hj.transpose() * S.inverse();
  
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
