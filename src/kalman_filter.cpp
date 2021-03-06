#include "kalman_filter.h"
#include "tools.h"
#include<math.h>
#include<iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

// 预测状态值
void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

/**
 *计算卡尔曼滤波
 **/
void KalmanFilter::Update(const VectorXd &z) {
  // 使用转换矩阵计算K值，K值可以判断预测值和测量值在最终结果的比例
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd St = S.inverse();
  MatrixXd K = P_ * Ht * St;
  // 更新状态和协方差矩阵
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

/**
 *计算扩展卡尔曼滤波
 **/
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // 使用雅各比矩阵计算K值，K值可以判断预测值和测量值在最终结果的比例
  Tools tools;
  MatrixXd Hj = tools.CalculateJacobian(x_);
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R_;
  MatrixXd St = S.inverse();
  MatrixXd K = P_ * Hjt * St;
  // 把笛卡尔坐标转换成极坐标，转换成和测量值相同的格式
  VectorXd z_pred = VectorXd(3);
  z_pred[0] = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
  z_pred[1] = atan2(x_[1], x_[0]);
  if (z_pred[0] != 0) {
    z_pred[2] = (x_[0] * x_[2] + x_[1] * x_[3]) / z_pred[0];
  } else {
    z_pred[2] = 0;
  }
  // 更新状态和协方差矩阵
  VectorXd y = z - z_pred;
  // 注意：这个地方需要把角度范围缩放到-PI到PI之间
  while (y[1] > M_PI){
    y[1] -= 2 * M_PI;
  }
  while (y[1] < -M_PI){
    y[1] += 2 * M_PI;
  }
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
