#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
 is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  std::cout<<"Process Measurements " << std::endl;
  std::cout<<"Sensor Type"<<meas_package.sensor_type_<<std::endl;
  if (!is_initialized_) {
    //cout << "Kalman Filter Initialization " << endl;
    std::cout<<"Init " << std::endl;

    // set the state with the initial location and zero velocity
    if(use_radar_ && (meas_package.sensor_type_ == MeasurementPackage::RADAR))
    {
      std::cout<<"Init radar " << std::endl;
    x_ << meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]), 
          meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]), 
          meas_package.raw_measurements_[2],
      	  0,
          0;
     P_ << 1, 0, 0, 0, 0,
           0, 1, 0, 0, 0,
           0, 0, 1, 0, 0, 
           0, 0, 0, 1, 0, 
           0, 0, 0, 0, 1;
    
    }
    else if(use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)){
      std::cout<<"Init laser " << std::endl;
    x_ << meas_package.raw_measurements_[0], 
          meas_package.raw_measurements_[1], 
          0,
      	  0,
          0;
    
    P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
      	  0, std_laspy_*std_laspy_, 0, 0, 0,
          0, 0, 5, 0, 0, 
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    }

    previous_timestamp_  = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
   
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);
  if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  { UpdateLidar(meas_package);}
  else if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {UpdateRadar(meas_package);}
  return;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
   // set state dimension
 
   std::cout<<"Prediction" << std::endl;
 
   n_x_ = 5;

  // set augmented dimension
  n_aug_ = 7;

  // Process noise standard deviation longitudinal acceleration in m/s^2
   std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
   std_yawdd_ = 1.0;

  // define spreading parameter
   lambda_ = 3 - n_aug_;
  
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  /**
   * Student part begin
   */
 
  // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug[n_x_ ] = 0;
  x_aug[n_x_ + 1] = 0;

  // create augmented covariance matrix
  P_aug.fill(0);
  P_aug.topLeftCorner(P_.rows(),P_.cols()) = P_;
  P_aug(P_.rows(),P_.cols()) = std_a_*std_a_;
  P_aug(P_.rows()+1,P_.cols()+1) = std_yawdd_*std_yawdd_;
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  for(int i = 0;i< n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i);
    Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
  }

  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    // write predicted sigma points into right column
  
  for(int sig_n = 0; sig_n < (2*n_aug_ + 1); sig_n++)
  {
   if(fabs(Xsig_aug(4,sig_n)) > 0.001) { 
        Xsig_pred_(0,sig_n) = Xsig_aug(0,sig_n) + (Xsig_aug(2,sig_n)/Xsig_aug(4,sig_n))*(sin(Xsig_aug(3,sig_n)+Xsig_aug(4,sig_n)*delta_t)-sin(Xsig_aug(3,sig_n))) + (0.5*delta_t*delta_t*cos(Xsig_aug(3,sig_n)*Xsig_aug(5,sig_n)));
        Xsig_pred_(1,sig_n) = Xsig_aug(1,sig_n) + (Xsig_aug(2,sig_n)/Xsig_aug(4,sig_n))*(-cos(Xsig_aug(3,sig_n)+Xsig_aug(4,sig_n)*delta_t)+cos(Xsig_aug(3,sig_n))) + (0.5*delta_t*delta_t*sin(Xsig_aug(3,sig_n)*Xsig_aug(5,sig_n)));
   }
   else
      {
        Xsig_pred_(0,sig_n) = Xsig_aug(0,sig_n) + (Xsig_aug(2,sig_n)*delta_t*cos(Xsig_aug(2,sig_n))) +(0.5*delta_t*delta_t*cos(Xsig_aug(3,sig_n)*Xsig_aug(5,sig_n))) ;
        Xsig_pred_(1,sig_n) = Xsig_aug(1,sig_n) + (Xsig_aug(2,sig_n)*delta_t*sin(Xsig_aug(2,sig_n))) + (0.5*delta_t*delta_t*sin(Xsig_aug(3,sig_n)*Xsig_aug(5,sig_n)));
      }
        Xsig_pred_(2,sig_n) = Xsig_aug(2,sig_n) + delta_t*Xsig_aug(5,sig_n);
        Xsig_pred_(3,sig_n) = Xsig_aug(3,sig_n) +  Xsig_aug(4,sig_n)*delta_t +0.5*delta_t*delta_t*Xsig_aug(6,sig_n);
        Xsig_pred_(4,sig_n) = Xsig_aug(4,sig_n) + delta_t*Xsig_aug(6,sig_n);
      
  }
  
   // create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  
  // set weights
  weights(0) = lambda_/(lambda_ +n_aug_);
  
  for(int i=1; i<(2*n_aug_+1);i++)
    weights(i) = 1/(2*(lambda_ +n_aug_));

  // predict state mean
  x_.fill(0.0);
  for(int state_n = 0; state_n <n_x_; state_n++)
  {
      
      for(int sig_n =0; sig_n < (2*n_aug_+1); sig_n++)
      {
          x_(state_n) = x_(state_n) + weights(sig_n)* Xsig_pred_(state_n,sig_n);
      }
  
  }
  // predict state covariance matrix
  
    P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
  }
  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
   std::cout<<"Lidar update" << std::endl;
 
    VectorXd z{ meas_package.raw_measurements_ };
  
   MatrixXd H_ { MatrixXd(2, n_x_) };
    H_  << 1, 0, 0, 0, 0,
           0, 1, 0, 0, 0;
  
   MatrixXd R_{ MatrixXd(2, 2) };
   R_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;
  
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // set measurement dimension, radar can measure r, phi, and r_dot
  std::cout<<"Radar Update" << std::endl;

  int n_z_ = 3;
  
  // set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  double weight = 0.5/(lambda_+n_aug_);
  weights(0) = weight_0;

  for (int i=1; i<2*n_aug_+1; ++i) {  
    weights(i) = weight;
  }
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  
  // transform sigma points into measurement space
  for(int sig_n = 0; sig_n < (2*n_aug_+1); sig_n++)
  {
      Zsig(0,sig_n) = sqrt((Xsig_pred_(0,sig_n)*Xsig_pred_(0,sig_n))+(Xsig_pred_(1,sig_n)*Xsig_pred_(1,sig_n)));
      Zsig(1,sig_n) = atan2(Xsig_pred_(1,sig_n),Xsig_pred_(0,sig_n));
      Zsig(2,sig_n) = ((Xsig_pred_(0,sig_n)*cos(Xsig_pred_(3,sig_n))*Xsig_pred_(2,sig_n))+(Xsig_pred_(1,sig_n)*sin(Xsig_pred_(3,sig_n))*Xsig_pred_(2,sig_n)))/Zsig(0,sig_n);
  }
  
  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }


  // calculate innovation covariance matrix S
      S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
      // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    S = S + weights(i) * z_diff * z_diff.transpose() ;
  }
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R<< std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
  S=S+R;
  
   VectorXd z{ meas_package.raw_measurements_ };
  
   // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  
  Tc.fill(0.0);

  for(int sig_n = 0; sig_n < 2*n_aug_ + 1; sig_n++)
  {
       
       MatrixXd x_diff = Xsig_pred_.col(sig_n) - x_;
     
       MatrixXd z_diff = Zsig.col(sig_n) - z_pred;
       
       while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
       while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    
 
       // angle normalization
       while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
       while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
       
       Tc = Tc + weights(sig_n)*x_diff*z_diff.transpose();
       
  }

  // calculate Kalman gain K;
  
  MatrixXd K = Tc*S.inverse();

  // update state mean and covariance matrix
  VectorXd z_angle_difference{ z - z_pred };
  
  while (z_angle_difference(1)> M_PI) z_angle_difference(1)-=2.*M_PI;
  while (z_angle_difference(1)<-M_PI) z_angle_difference(1)+=2.*M_PI;
  
  x_ = x_ + K*(z_angle_difference);
  P_ = P_ - K*S*K.transpose();
  
}