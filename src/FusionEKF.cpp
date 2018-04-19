#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // 初始化一些计算中需要用到的矩阵
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  P_ = MatrixXd(4, 4);
  F_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);

  // 激光雷达测量误差协方差
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // 雷达测量误差协方差
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // 初始化激光雷达的转换矩阵
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // 初始化雷达的雅各比矩阵
  Hj_ << 1, 0, 0, 0,
         0, 1, 0, 0,
		     0, 0, 1, 0;

  // 初始化P矩阵
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
  
  // 初始化F矩阵
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
  
  // 初始化Q矩阵
  Q_ <<  0, 0, 0, 0, 
		     0, 0, 0, 0,
		     0, 0, 0, 0,
		     0, 0, 0, 0;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  初始化，使用第一次的测量数据去更新状态且初始化需要用到的协方差矩阵
   ****************************************************************************/
  if (!is_initialized_) {
    VectorXd x_ = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // 雷达测量数据是极坐标数据，转换为笛卡尔坐标，这样我们的模型就和激光雷达模型一致
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];
      x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);
      // 各个参数初始化
      ekf_.Init(x_, P_, F_, Hj_, R_radar_, Q_);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // 激光雷达测量数据只有位移，没有速度，初始化速度为0
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      // 各个参数初始化
      ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
    }
    // 初始化时间
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }


  /*****************************************************************************
   *  预测，根据时间差和加速度去预测物体的位置和速度
   ****************************************************************************/
  // 计算两次测量时间差
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  // 在F添加时间融合
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  float noise_ax = 9;
  float noise_ay = 9;
  // 修改预测误差协方差
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0, 
             0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
             dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
             0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;
  ekf_.Predict();
  
  
  /*****************************************************************************
   *  更新，使用测量数据更新状态和P矩阵
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
