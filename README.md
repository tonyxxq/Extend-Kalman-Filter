- 项目简介

  执行这个项目结果需要先下载term2的模拟器，[地址](https://github.com/udacity/self-driving-car-sim/releases)

  因为需要使用uWebSocketIO和模拟器进行连接，在项目下有[install-mac.sh](./install-mac.sh)和[install-ubuntu.sh](./install-ubuntu.sh)两个文件，分别针对不同的平台安装uWebSocketIO。

  本项目在mac上使用eclipse开发。


    输入: 模拟器产生的数据，包括激光雷达和雷达数据。

    输出: c++代码输出数据，该数据传输到模拟器。

- 编译环境要求

  - cmake >= 3.5
  - make >= 4.1 (Linux, Mac), 3.81 (Windows)
  - gcc/g++ >= 5.4

- 程序下载,编译,执行过程

  Clone this repo.

  `mkdir build && cd build`

  `cmake .. && make` 

  `./ExtendedKF `

- 部分代码

  > 执行过程
  >
  > 1. 初始化数据
  >    - 激光雷达直接使用第一次的测量数据更新就可以了。
  >    - 毫米波雷达需要把第一次的测量数据从极坐标转换成笛卡尔坐标的数据，这样激光雷达和毫米波雷达的数据格式一致，使用同一个预测函数。
  > 2. 预测，
  >    - 激光雷达和毫米波雷达使用同一个预测函数，预测的时候一共有两个变量值需要更新到F合Q矩阵，时间段和加速度。
  > 3. 更新
  >    - 激光雷达使用普通的卡尔曼滤波更新过程
  >    - 因为卡尔曼滤波只能满足高斯分布，但是非线性函数执行后，结果就不是线性的了。毫米波雷达的测量值不满足线性条件，需要使用扩展卡尔曼滤波，把测量函数转换为近似的线性函数，再处理。

  - FusionEKF.cpp

    ```c++
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
    ```

  - kalman_filter.cpp

    ```c++
    #include "kalman_filter.h"
    #include "tools.h"

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
      MatrixXd Si = S.inverse();
      MatrixXd K = P_ * Ht* Si;
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
      MatrixXd Si = S.inverse();
      MatrixXd K = P_ * Hjt * Si;
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
      x_ = x_ + K * y;
      long x_size = x_.size();
      MatrixXd I = MatrixXd::Identity(x_size, x_size);
      P_ = (I - K * Hj) * P_;
    }

    ```

  - tools.cpp

    ```c++
    #include <iostream>
    #include "tools.h"

    using Eigen::VectorXd;
    using Eigen::MatrixXd;
    using std::vector;

    Tools::Tools() {
    }

    Tools::~Tools() {
    }

    /**
     * 计算RMS，测试预测值和真实值的误差
     */
    VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                           const vector<VectorXd> &ground_truth) {
      VectorXd rmse(4);
      rmse << 0, 0, 0, 0;

      // check the validity of the following inputs:
      //  * the estimation vector size should not be zero
      //  * the estimation vector size should equal ground truth vector size
      if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
      }

      //accumulate squared residuals
      for (unsigned int i = 0; i < estimations.size(); ++i) {

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array() * residual.array();
        rmse += residual;
      }

      //calculate the mean
      rmse = rmse / estimations.size();

      //calculate the squared root
      rmse = rmse.array().sqrt();

      //return the result
      return rmse;
    }

    /**
     * 计算雅各比矩阵
     */
    MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

      MatrixXd Hj(3, 4);
      //recover state parameters
      float px = x_state(0);
      float py = x_state(1);
      float vx = x_state(2);
      float vy = x_state(3);

      //pre-compute a set of terms to avoid repeated calculation
      float c1 = px * px + py * py;
      float c2 = sqrt(c1);
      float c3 = (c1 * c2);

      //check division by zero
      if (fabs(c1) < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
      }

      //compute the Jacobian matrix
      Hj << (px / c2), (py / c2), 0, 0, -(py / c1), (px / c1), 0, 0, py
          * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py
          / c2;

      return Hj;
    }
    ```

    ​


