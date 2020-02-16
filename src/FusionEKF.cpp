#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  ekf_.P_ = MatrixXd::Zero(4,4);
  ekf_.P_(0,0) = ekf_.P_(1,1) = 1;
  ekf_.P_(2,2) = ekf_.P_(3,3) = 1000;
  
  H_laser_ << 1,0,0,0,
              0,1,0,0;
  
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    std::cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      ekf_.x_[0] = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
      ekf_.x_[1] = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
      ekf_.x_[2] = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
      ekf_.x_[3] = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);
      if(ekf_.x_[0]< 0.0001) ekf_.x_[0] = 0.0001;
      if(ekf_.x_[1]< 0.0001) ekf_.x_[1] = 0.0001;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_[0] = measurement_pack.raw_measurements_[0];
      ekf_.x_[1] = measurement_pack.raw_measurements_[1];
      ekf_.x_[2] = 0;
      ekf_.x_[3] = 0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_ ;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  
  MatrixXd G(4,2);
  G(0,1) = G(1,0) = G(2,1) = G(3,0) = 0;
  G(1,1) = G(0,0) = dt*dt/2.0;
  G(2,0) = G(3,1) = dt;
  
  MatrixXd Q(2,2);
  Q(0,1) = Q(1,0) = 0;
  Q(0,0) = Q(1,1) = 9;
  ekf_.Q_ = G*Q*G.transpose();
  
  ekf_.F_ = MatrixXd::Identity(4, 4);
  ekf_.F_(0,2) = ekf_.F_(1,3) = dt;
  
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  previous_timestamp_ = measurement_pack.timestamp_;
  // print the output
//   std::cout << "x_ = " << ekf_.x_ << std::endl;
//   std::cout << "P_ = " << ekf_.P_ << std::endl;
}
