#include "FusionEKF.h"
#include "tools.h"
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.F_ = MatrixXd::Zero(4,4);
  ekf_.F_(0,0) = 1.0;
  ekf_.F_(1,1) = 1.0;
  ekf_.F_(2,2) = 1.0;
  ekf_.F_(3,3) = 1.0;
  
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;
              
  ekf_.H_ = MatrixXd(2,4);
  ekf_.H_ << 1,0,0,0,
             0,1,0,0;
             
  ekf_.Q_ = MatrixXd::Zero(4,4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    ekf_.Q_ = MatrixXd::Zero(4,4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd raMeas = measurement_pack.raw_measurements_;
      if ((raMeas(0) < 0.001) && (raMeas(1) < 0.001) && (raMeas(2) < 0.001))
        return;
        
      ekf_.x_[0] = raMeas[0] * cos(raMeas[1]);
      ekf_.x_[1] = raMeas[0] * sin(raMeas[1]);
      ekf_.x_[2] = raMeas[2] * cos(raMeas[1]);
      ekf_.x_[3] = raMeas[2] * sin(raMeas[1]);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      VectorXd laMeas = measurement_pack.raw_measurements_;
      if ((laMeas(0) < 0.001) && (laMeas(1)  < 0.001))
          return;
      ekf_.x_ << laMeas(0), laMeas(1), 0, 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
   
    /*
    cout<<"x_:\n" << ekf_.x_<<endl;
    cout<<"P_:\n" << ekf_.P_<<endl;
    cout<<"F_:\n" << ekf_.F_<<endl;
    cout<<"H_:\n" << ekf_.H_<<endl;
    cout<<"Q_:\n" << ekf_.Q_<<endl;
    
    cout<<"\n"<<endl;
    */
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  double delta = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  //cout<<"delta:"<<delta<<endl;

  /* if delta is zero, skip predict */
  if (delta > 1e-6) {
    ekf_.F_(0,2) = delta;
    ekf_.F_(1,3) = delta;
  
    //double noise_ax2 = 9*9;
    double noise_ax2 = 9;
    //double noise_ay2 = 9*9;
    double noise_ay2 = 9;
    
    ekf_.Q_(0,0) = pow(delta, 4)/4 * noise_ax2;
    ekf_.Q_(0,2) = pow(delta, 3)/2 * noise_ax2;
    ekf_.Q_(1,1) = pow(delta, 4)/4 * noise_ay2;
    ekf_.Q_(1,3) = pow(delta, 3)/2 * noise_ay2;
    ekf_.Q_(2,0) = pow(delta, 3)/2 * noise_ax2;
    ekf_.Q_(2,2) = pow(delta, 2) * noise_ax2;
    ekf_.Q_(3,1) = pow(delta, 3)/2 * noise_ay2;
    ekf_.Q_(3,3) = pow(delta, 2) * noise_ay2;
 
    ekf_.Predict();
  }
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  
  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << "\n" <<endl;
  
}