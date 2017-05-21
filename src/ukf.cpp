#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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

  P_ << 1,0,0,0,0,
          0,1,0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* State dimension
  n_x_ = x_.size();

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;

  cnt = 0;

  //create example vector for mean predicted measurement
  z_pred_ = VectorXd(n_z_);


  S_ = MatrixXd(n_z_, n_z_);

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  // number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // predicted sigma points for x state
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  ///* Weights of sigma points
  weights_ = VectorXd(n_sig_);

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* the current NIS for radar
  NIS_radar_ = 0;

  ///* the current NIS for laser
  NIS_laser_ = 0;

  //measurement noise initialization
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_ , 0, 0,
             0, std_radphi_*std_radphi_,0,
             0, 0, std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
             0, std_laspy_*std_laspy_;

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

  if (!is_initialized_){

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      float px = rho*cos(phi);
      float py = rho*sin(phi);
      float vx = rho_dot*cos(phi);
      float vy = rho_dot*sin(phi);
      float v = sqrt(vx*vx + vy*vy);

      x_ << px, py, v, 0, 0;


    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

    }

  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); i++){
    weights_(i) = 0.5 / (lambda_ + n_aug_);
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  double dt_ = meas_package.timestamp_ - time_us_;

  // time unit to seconds
  dt_ /= 1000000.0;

  time_us_ = meas_package.timestamp_;

  Prediction(dt_);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
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

  VectorXd x_aug = VectorXd(n_aug_);

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  x_aug.fill(0.0);

  x_aug.head(n_x_) = x_;

  P_aug.fill(0);

  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  P_aug(5,5) = std_a_*std_a_;

  P_aug(6,6) = std_yawdd_*std_yawdd_;


  MatrixXd L = P_aug.llt().matrixL();


  //  Generated sigma points with Xsig_aug (px, py, v, yaw, yawd, nu_a, nu_yawdd)
  double sqrt_lambda_n_aug = sqrt(lambda_+n_aug_); // Save some computations
  VectorXd sqrt_lambda_n_aug_L;

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++) {
    sqrt_lambda_n_aug_L = sqrt_lambda_n_aug * L.col(i);
    Xsig_aug.col(i+1)        = x_aug + sqrt_lambda_n_aug_L;
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda_n_aug_L;
  }


  // Predict sigma points with Xsig_aug (px, py, v, yaw, yawd) with generated sigma points(px, py, v, yaw, yawd, nu_a, nu_yawdd)
  for (int i = 0; i<n_sig_; i++){
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double p_xp, p_yp;

    if (fabs(yawd) > 0.001) {
      p_xp = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      p_yp = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      p_xp = p_x + v*delta_t*cos(yaw);
      p_yp = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

//    double v_p = v + delta_t*nu_a;
//    double yaw_p =  + 0.5*delta_t*delta_t*nu_yawdd;
//    double yawd_p = yawd + delta_t*nu_yawdd;


    p_xp += 0.5*(delta_t*delta_t)*cos(yaw)*nu_a;
    p_yp += 0.5*(delta_t*delta_t)*sin(yaw)*nu_a;
    v_p += nu_a*delta_t;
    yaw_p += 0.5*delta_t*delta_t*nu_yawdd;
    yawd_p += delta_t*nu_yawdd;

    Xsig_pred_(0, i) = p_xp;
    Xsig_pred_(1, i) = p_yp;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;

  }

  // Predicted state mean [5x1] - vectorized sum
  x_ = Xsig_pred_ * weights_;


  // Predicted Covariance Matrix [5x1] iterate over sigma points
  P_.fill(0.0);

  for (int i = 0; i < n_sig_; i++){

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

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

  int n_z = 2;

  // Create matrix for sigma points in measurement space
  // Transform sigma points into measurement space

  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);
  UpdateUKF(meas_package, Zsig, n_z);

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

  double n_z = 3;

  MatrixXd Zsig = Xsig_pred_.block(0,0,n_z, n_sig_);

  for (int i = 0; i<n_sig_; i++){

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;


    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
//    if (p_x<0.00001){
//        Zsig(1,i) = atan2(p_y,0.00001);
//    }
//    else{
    Zsig(1,i) = atan2(p_y,p_x);
//    }
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);

  }

  UpdateUKF(meas_package, Zsig, n_z);

}


void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z){
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred  = Zsig * weights_;
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    R = R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    R = R_lidar_;
  }
  S = S + R;

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
      // Angle normalization
      NormAng(&(z_diff(1)));
    }

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    NormAng(&(x_diff(3)));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // Measurements
  VectorXd z = meas_package.raw_measurements_;
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // Residual
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    // Angle normalization
    NormAng(&(z_diff(1)));
  }
  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  // Calculate NIS


  cout<<"total count "<<cnt<<endl;

  cnt += 1;


  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    NIS_radar_ = z.transpose() * S.inverse() * z;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    NIS_laser_ = z.transpose() * S.inverse() * z;
  }

}

void UKF::NormAng(double *ang) {
    while (*ang > M_PI) *ang -= 2. * M_PI;
    while (*ang < -M_PI) *ang += 2. * M_PI;
}
