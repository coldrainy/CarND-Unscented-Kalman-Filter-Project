#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);

  weights_ = VectorXd(2*n_aug_+1);

  use_laser_ = true;

  use_radar_ = true;

  is_initialized_ = false;
}
UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if(is_initialized_ == false){
        x_.fill(0);
        P_<<0.15,0,0,0,0,
           0,0.15,0,0,0,
           0,0,1,0,0,
           0,0,0,1,0,
           0,0,0,0,1;

        if(meas_package.sensor_type_ == MeasurementPackage::LASER){
            float x = meas_package.raw_measurements_(0);
            float y = meas_package.raw_measurements_(1);
            x_(0) = x;
            x_(1) = y;
        }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            float ro = meas_package.raw_measurements_(0);
            float theta = meas_package.raw_measurements_(1);
            float x = ro*cos(theta);
            float y = ro*sin(theta);
            x_(0) = x;
            x_(1) = y;
        }
        std::cout<<"------------initialize-------------"<<std::endl;
        std::cout<<"measurement:"<<std::endl;
        std::cout<<meas_package.raw_measurements_<<std::endl;
        std::cout<<"x_:"<<std::endl;
        std::cout<<x_<<std::endl;
        std::cout<<"P_:"<<std::endl;
        std::cout<<P_<<std::endl;


        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }

    long long current_time = meas_package.timestamp_;
    float dt = (current_time - time_us_)/1000000.0;
    time_us_ = current_time;


    Prediction(dt);
    std::cout<<"------------predict-------------"<<std::endl;
    std::cout<<"dt:"<<dt<<std::endl;
    std::cout<<"x_:"<<std::endl;
    std::cout<<x_<<std::endl;
    std::cout<<"P_:"<<std::endl;
    std::cout<<P_<<std::endl;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true){
        UpdateLidar(meas_package);
        std::cout<<"------------updatelidar-------------"<<std::endl;
        std::cout<<"measurement:"<<std::endl;
        std::cout<<meas_package.raw_measurements_<<std::endl;

        std::cout<<"x_:"<<std::endl;
        std::cout<<x_<<std::endl;
        std::cout<<"P_:"<<std::endl;
        std::cout<<P_<<std::endl;

    }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_== true){
        UpdateRadar(meas_package);
        std::cout<<"------------updateradar-------------"<<std::endl;
        std::cout<<"measurement:"<<std::endl;
        std::cout<<meas_package.raw_measurements_<<std::endl;

        std::cout<<"x_:"<<std::endl;
        std::cout<<x_<<std::endl;
        std::cout<<"P_:"<<std::endl;
        std::cout<<P_<<std::endl;

    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

    VectorXd x_aug(n_aug_);
    x_aug.fill(0);
    x_aug.head(n_x_) = x_;


    MatrixXd P_aug(n_aug_,n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(n_x_,n_x_) = std_a_*std_a_;
    P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;


    MatrixXd Xsig_aug(n_aug_,2*n_aug_+1);
    Xsig_aug.col(0) = x_aug;
    MatrixXd A = P_aug.llt().matrixL();
    for(int i=0;i<n_aug_;i++){
        Xsig_aug.col(i+1) = x_aug+sqrt(lambda_+n_aug_)*A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug-sqrt(lambda_+n_aug_)*A.col(i);
    }
    for(int i=0;i<2*n_aug_+1;i++){
        float px = Xsig_aug(0,i);
        float py = Xsig_aug(1,i);
        float v = Xsig_aug(2,i);
        float yaw = Xsig_aug(3,i);
        float yaw_rate = Xsig_aug(4,i);
        float a = Xsig_aug(5,i);
        float yaw_acc = Xsig_aug(6,i);

        double px_p,py_p;
        if(fabs(yaw_rate)>0.001){
          px_p = px+v/yaw_rate*(sin(yaw+yaw_rate*delta_t)-sin(yaw));
          py_p = py+v/yaw_rate*(-cos(yaw+yaw_rate*delta_t)+cos(yaw));
        }else{
          px_p = px+v*cos(yaw)*delta_t;
          py_p = py+v*sin(yaw)*delta_t;
        }
        double v_p = v;
        double yaw_p = yaw+yaw_rate*delta_t;
        double yaw_rate_p = yaw_rate;

        px_p = px_p + 0.5*delta_t*delta_t*cos(yaw)*a;
        py_p = py_p + 0.5*delta_t*delta_t*sin(yaw)*a;
        v_p = v_p + delta_t*a;
        yaw_p = yaw_p + 0.5*delta_t*delta_t*yaw_acc;
        yaw_rate_p = yaw_rate_p+delta_t*yaw_acc;

        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yaw_rate_p;
    }
    for(int i=0;i<2*n_aug_+1;i++){
        if(i==0){
            double weight = lambda_/(lambda_+n_aug_);
            weights_(i) = weight;
        }
        else{
            double weight =  0.5/(lambda_+n_aug_);
            weights_(i) = weight;
        }
    }

    x_ = Xsig_pred_*weights_;

    P_.fill(0);
    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd diff = Xsig_pred_.col(i)-x_;
        while (diff(3)> M_PI) diff(3)-=2*M_PI;
        while (diff(3)<-M_PI) diff(3)+=2*M_PI;
        P_+=weights_(i)*diff*diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    VectorXd z_pred(2);
    MatrixXd Zsig(2,2*n_aug_+1);
    Zsig.fill(0);
    for(int i=0;i<2*n_aug_+1;i++){
        float px = Xsig_pred_(0,i);
        float py = Xsig_pred_(1,i);
        Zsig(0,i) = px;
        Zsig(1,i) = py;
    }
    z_pred =  Zsig * weights_;
    MatrixXd S(2,2);
    S.fill(0);
    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * diff * diff.transpose();
    }
    MatrixXd R(2,2);
    R << std_laspx_*std_laspx_,0,
         0,std_laspy_*std_laspy_;
    S+=R;

    MatrixXd T(n_x_,2);
    T.fill(0);
    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        if(x_diff(3)>M_PI) x_diff(3) -=2* M_PI;
        if(x_diff(3)<-M_PI) x_diff(3) +=2* M_PI;
        VectorXd z_diff = Zsig.col(i) - z_pred;
        T = T+weights_(i)*x_diff*z_diff.transpose();
    }

    MatrixXd K(n_x_,2);
    K = T*S.inverse();

    VectorXd z(2);
    z<<meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);
    VectorXd z_diff = z-z_pred;
    x_ = x_+K*z_diff;
    P_ = P_ - K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    VectorXd z_pred(3);
    MatrixXd Zsig(3,2*n_aug_+1);
    Zsig.fill(0);

    for(int i=0;i<2*n_aug_+1;i++){
        float px = Xsig_pred_(0,i);
        float py = Xsig_pred_(1,i);
        float v = Xsig_pred_(2,i);
        float phi = Xsig_pred_(3,i);

        Zsig(0,i) = sqrt(px*px+py*py);
        if(px < 0.001 && py < 0.001){
            std::cout<<" line 325 error"<<std::endl;
            return;
        }
        Zsig(1,i) = atan2(py,px);
        Zsig(2,i) = (px*v*cos(phi)+py*v*sin(phi))/Zsig(0,i);
    }
    z_pred = Zsig*weights_;

    MatrixXd S(3,3);
    S.fill(0);
    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        if(z_diff(1)>M_PI) z_diff(1)-=2*M_PI;
        if(z_diff(1)<-M_PI) z_diff(1)+=2*M_PI;
        S += weights_(i)*z_diff*z_diff.transpose();
    }
    MatrixXd R(3,3);
    R<<std_radr_*std_radr_,0,0,
        0,std_radphi_*std_radphi_,0,
        0,0,std_radrd_*std_radrd_;
    S += R;
    MatrixXd T(n_x_,3);
    T.fill(0);
    for(int i = 0;i<2*n_aug_+1;i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        if(x_diff(3)>M_PI) x_diff(3) -= 2*M_PI;
        if(x_diff(3)<-M_PI) x_diff(3) += 2*M_PI;
        VectorXd z_diff = Zsig.col(i) - z_pred;
        if(z_diff(1)>M_PI) z_diff(1)-=2*M_PI;
        if(z_diff(1)<-M_PI) z_diff(1)+=2*M_PI;
        T += weights_(i)*x_diff*z_diff.transpose();
    }
    MatrixXd K(n_x_,3);
    K = T*S.inverse();

    VectorXd z(3);
    z<<meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);
    VectorXd z_diff = z-z_pred;
    if(z_diff(1)>M_PI) z_diff(1)-=2*M_PI;
    if(z_diff(1)<-M_PI) z_diff(1)+=2*M_PI;

    x_ = x_ + K*z_diff;
    P_ = P_ - K*S*K.transpose();
}
