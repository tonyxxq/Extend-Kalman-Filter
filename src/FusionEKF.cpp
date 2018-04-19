#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
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
  //state covariance matrix P
  ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

  //the initial transition matrix F_
  ekf_.F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;
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
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // 测量数据是极坐标数据，转换为笛卡尔坐标
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // 雷达测量数据只有位移，没有速度，初始化速度为0
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    // 初始化时间
    previous_timestamp_ = measurement_pack.timestamp_;
    // 设置是否已初始化为true
    is_initialized_ = true;
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
  // 计算时间差
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float noise_ax = 9;
  float noise_ay = 9;
  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0, 0, dt_4 / 4
      * noise_ay, 0, dt_3 / 2 * noise_ay, dt_3 / 2 * noise_ax, 0, dt_2
      * noise_ax, 0, 0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    // Radar updates
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
