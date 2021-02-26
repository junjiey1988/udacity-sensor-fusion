#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

UKF::UKF() {
  time_us_ = 0;

  is_initialized_ = false;


  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI;

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
  
  // Sigma points parameters
  n_aug_ = 7;
  n_sig_ = 2 * n_aug_ + 1; 
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  lambda_ = 3 - n_aug_;

  // Sigma points weights
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  long dt_us = meas_package.timestamp_ - time_us_;
  time_us_ = meas_package.timestamp_;

  if (is_initialized_)
  {
    Prediction(dt_us / 1000000.0);
  }

  // Updates Lidar measurement.
  if (meas_package.sensor_type_ == meas_package.LASER)
  {
    if (!is_initialized_)
    {
      x_.fill(0.0);
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

      is_initialized_ = true;
    }
    else if (use_laser_)
    {
      UpdateLidar(meas_package);
    }
  }
  // Updates Radar measurement.
  else if (meas_package.sensor_type_ == meas_package.RADAR && use_radar_ && is_initialized_)
  {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  // Generates sigma points 
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL(); // square root of P
  

  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // Predicts sigma points 
  for (int i = 0; i < n_sig_; ++i)
  {
    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    if (fabs(yawd) > std::numeric_limits<double>::epsilon()) 
    {
      // Prevents division by zero
      px = px + (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
      py = py + (v / yawd) * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      // Uses linear model when no rotation.
      px = px + (v * delta_t * cos(yaw));
      py = py + (v * delta_t * sin(yaw));
    }

    yaw = yaw + (yawd * delta_t);

    // Adds noise
    Xsig_pred_(0, i) = px + (0.5 * nu_a * delta_t * delta_t * cos(yaw));
    Xsig_pred_(1, i) = py + (0.5 * nu_a * delta_t * delta_t * sin(yaw));
    Xsig_pred_(2, i) = v + (nu_a * delta_t);
    Xsig_pred_(3, i) = yaw + (0.5 * nu_yawdd * delta_t * delta_t);
    Xsig_pred_(4, i) = yawd + (nu_yawdd * delta_t);
  }

  // Predicts mean and covariance
  x_.fill(0.0);
  P_.fill(0.0);

  for (int i = 0; i < n_sig_; ++i)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  for (int i = 0; i < n_sig_; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Normalizes angle
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.0 * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Predicts measurement 
  int n_z = 2;                      // measurement dimension
  MatrixXd H = MatrixXd(n_z, n_x_); // measurement matrix
  
  H.fill(0.0);
  H(0, 0) = 1;
  H(1, 1) = 1;

  VectorXd z_pred = H * x_;

  // Updates state 
  VectorXd z = VectorXd(n_z); // measurement
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);

  MatrixXd R = MatrixXd(n_z, n_z); // measurement covariance
  R.fill(0.0);
  R(0, 0) = std_laspx_ * std_laspx_;
  R(1, 1) = std_laspy_ * std_laspy_;

  VectorXd y = z - z_pred;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);

  x_ = x_ + (K * y);
  P_ = (I - K * H) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Predict measurement 
  int n_z = 3;                            // measurement dimension
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);  // measurement space sigma points
  VectorXd z_pred = VectorXd(n_z);        // mean predicted measurement
  MatrixXd S = MatrixXd(n_z, n_z);        // measurement covariance matrix

  for (int i = 0; i < n_sig_; i++)
  {
    double px  = Xsig_pred_(0, i);
    double py  = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    Zsig(0, i) = sqrt(px * px + py * py); 
    Zsig(1, i) = atan2(py, px); 
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig(0, i); 
  }

  z_pred.fill(0.0);
  S.fill(0.0);

  for (int i = 0; i < n_sig_; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  for (int i = 0; i < n_sig_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.0 * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S(0, 0) += std_radr_ * std_radr_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_ * std_radrd_;

  // Update state 
  MatrixXd Tc = MatrixXd(n_x_, n_z); // cross correlation
  VectorXd z = VectorXd(n_z); // measurement

  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  Tc.fill(0.0);

  for (int i = 0; i < n_sig_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.0 * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.0 * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse(); 
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) >  M_PI) z_diff(1) -= 2.0 * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}