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
	//std_a_ = 30;
	std_a_ = 1.5;

	// Process noise standard deviation yaw acceleration in rad/s^2

	//std_yawdd_ = 30;
	std_yawdd_ = 0.57;

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

	is_initialized_ = false;
	time_us_ = 0;
	n_x_ = 5;
	n_aug_ = 7;
	n_z_ = 3;
	lambda_ = 3 - n_aug_;

	Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

	weights_ = VectorXd(2*n_aug_+1);
	double wt = lambda_/(lambda_+n_aug_);
	weights_(0) = wt;

	for (int i=1; i<2*n_aug_+1; i++) {  
		wt = 0.5/(n_aug_+lambda_);
		weights_(i) = wt;
	}

	R_rad_ = MatrixXd(3, 3);
	R_rad_ <<	std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0,std_radrd_*std_radrd_;

	R_lid_ = MatrixXd(2, 2);
	R_lid_ <<	std_laspx_*std_laspx_,0,
		0,std_laspy_*std_laspy_;
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

	if (!is_initialized_)
	{
		cout<<"Initilizing UKF...";

		x_ = VectorXd(n_x_);
		x_ << 0,0,0,0,0;

		P_ << 1,0,0,0,0,
			0,1,0,0,0,
			0,0,1,0,0,
			0,0,0,1,0,
			0,0,0,0,1;

		if (MeasurementPackage::LASER == meas_package.sensor_type_)
		{
			float px = meas_package.raw_measurements_(0);
			float py = meas_package.raw_measurements_(1);
			x_(0) = px;
			x_(1) = py;
			if (fabs(x_(0)) < 0.001 && fabs(x_(1)) < 0.001){
				x_(0) = 0.001;
				x_(1) = 0.001;
			}
			cout<<"with LASER data...";
		}
		else if (MeasurementPackage::RADAR == meas_package.sensor_type_)
		{
			float ro = meas_package.raw_measurements_(0);
			float phi = meas_package.raw_measurements_(1);
			float rod = meas_package.raw_measurements_(2);
			double vx = rod*cos(phi);
			double vy = rod*sin(phi);

			x_(0) = ro * cos(phi);
			x_(1) = ro * sin(phi);
			x_(2) = sqrt((vx*vx) + (vy*vy));
			cout<<"with RADAR data...";
		}
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		cout<<"done\n";
	}
	else 
	{
		double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
		time_us_ = meas_package.timestamp_;
		this->Prediction(dt);
		if(MeasurementPackage::LASER == meas_package.sensor_type_)
		{
			cout<<"LIDAR Update ";
			this->UpdateLidar(meas_package);
		}
		else
		{
			cout<<"RADAR Update ";
			this->UpdateRadar(meas_package);
		}
		//cout << "\nx_ =\t" << x_(0) << "\t" << x_(1) << "\t" << x_(2) << "\t" << x_(3) << "\t" << x_(4) << "\t" << endl;
		//cout << "P_ = \n" << P_ << endl;
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
	cout<<"Prediction...";
	//Create augmented state vector
	VectorXd x_aug_ = VectorXd(n_aug_);
	//Create augmented covariances matrix
	MatrixXd p_aug_ = MatrixXd(n_aug_, n_aug_);
	//Create augmented sigma point matrix
	MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	x_aug_.fill(0.0);
	x_aug_.head(5) = x_;
	x_aug_(5) = 0;
	x_aug_(6) = 0;

	p_aug_.fill(0.0);
	p_aug_.topLeftCorner(5,5) = P_;
	p_aug_(5,5) = std_a_*std_a_;
	p_aug_(6,6) = std_yawdd_*std_yawdd_;

	//Create square root matrix
	MatrixXd L = p_aug_.llt().matrixL();

	//Create augmented sigma points
	Xsig_aug_.col(0) = x_aug_;

	for(int i = 0; i < n_aug_; i++)
	{
		Xsig_aug_.col(i+1)       = x_aug_ + (sqrt(lambda_+n_aug_) * L.col(i));
		Xsig_aug_.col(i+1+n_aug_) = x_aug_ - (sqrt(lambda_+n_aug_) * L.col(i));  
	}
	cout<<"Sigma...";
	//Sigma point prediction
	for (int i = 0; i< 2*n_aug_+1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug_(0,i);
		double p_y = Xsig_aug_(1,i);
		double v = Xsig_aug_(2,i);
		double yaw = Xsig_aug_(3,i);
		double yawd = Xsig_aug_(4,i);
		double nu_a = Xsig_aug_(5,i);
		double nu_yawdd = Xsig_aug_(6,i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;
	}
	cout<<"State Mean...";
	//predicted state mean
	x_ = Xsig_pred_ * weights_; 
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		// State difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		// Angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
	}


	cout<<"Pred Done ";
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
	n_z_ = 2;
	//Predit measurement

	// Transform sigma points into measurement space
	MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z_, 2*this->n_aug_+1);

	VectorXd z_pred = VectorXd(n_z_);
	z_pred  = Zsig * weights_;
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);

	S.fill(0.0);
	for (int i = 0; i < 2*n_aug_+1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R_lid_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	Tc.fill(0.0);
	for (int i = 0; i < 2*n_aug_+1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}
	// Measurements
	VectorXd z = meas_package.raw_measurements_;
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
	NIS_laser_ = z.transpose() * S.inverse() * z;
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
	n_z_ = 3;
	//Predict Measurement
	MatrixXd Zsig = MatrixXd(this->n_z_, 2*this->n_aug_+1);

	//Tranform sigma points into measurement space
	for(int i=0; i<2*this->n_aug_+1;i++)
	{

		// extract values for better readibility
		double p_x = this->Xsig_pred_(0,i);
		double p_y = this->Xsig_pred_(1,i);
		double v  = this->Xsig_pred_(2,i);
		double yaw = this->Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig(1,i) = atan2(p_y,p_x);                                 //phi
		Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(this->n_z_);
	z_pred  = Zsig * weights_;

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(this->n_z_,this->n_z_);
	S.fill(0.0);
	for (int i = 0; i < 2 * this->n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + this->weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R_rad_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);
	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// Measurements
	VectorXd z = meas_package.raw_measurements_;

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	NIS_radar_ = z.transpose() * S.inverse() * z;
}
