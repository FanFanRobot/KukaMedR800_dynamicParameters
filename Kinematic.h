#pragma once
#include <math.h>
#include <iostream>
#include "mkl.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include<boost/progress.hpp>
#include<boost/timer.hpp>
#define _USE_MATH_DEFINES
using namespace Eigen;
using namespace std;
class Kinematic
{
public:
	Kinematic();
	~Kinematic();
	MatrixXd TransformOfAdiacentLink(double jointAngle, double alpha, double a, double d);
	MatrixXd Geometry_Jacobian_DH(VectorXd angles);
	MatrixXd PseudoInverse_JacobianGeometry(MatrixXd JacobianGeometry);
	void PseudoInverseSVD_JacobianGeometry(MatrixXd JacobianTrans, MatrixXd* M_pinv_);
	MatrixXd Geometry_JacobianTimeDerivative(VectorXd joints, VectorXd jointVelocity);
	MatrixXd Kinematic_KukaLBRMedR800(double angles[]);
	MatrixXd Kinematic_KukaLBRMedR800(VectorXd angles);
	VectorXd JointVelocity(VectorXd dj, VectorXd dj_);
	VectorXd CartesianVelocity(VectorXd vj, MatrixXd Jacobian_geometry);
	VectorXd OrientationError_matrix(MatrixXd desiredPose, MatrixXd currentPose);
	Vector3d OrientationError_quaternion(const Quaterniond &orientation_d, Quaterniond orientation);
	VectorXd PosiError(MatrixXd desiredPose, MatrixXd currentPose);
	VectorXd VelError(VectorXd desiredVel, VectorXd currentVel);
	double degree2Radian = M_PI / 180.f;
	//DH parameters of KUKA LBR Med R800
	double alpha0 = 0, alpha1 = -90.f*(degree2Radian), alpha2 = 90.f*(degree2Radian), alpha3 = 90.f*(degree2Radian),
		alpha4 = -90.f*(degree2Radian), alpha5 = -90.f*(degree2Radian), alpha6 = 90.f*(degree2Radian);
	double a = 0;
	double d1 = 0.34, d2 = 0, d3 = 0.4, d4 = 0, d5 = 0.4, d6 = 0, d7 = 0.126;
	Quaterniond orientation_d_{ Eigen::Quaterniond::Identity() };

	/*% Kinematics 	theta[rad] 	d[m]     a[m]      alpha[rad]
	% Joint 1       q1      0.34       0        	0
	% Joint 2       q2      0          0       - Π / 2
	% Joint 3       q3      0.4 	   0         Π / 2
	% Joint 4       q4      0          0 	     Π / 2
	% Joint 5       q5      0.4        0       - Π / 2
	% Joint 6       q6      0          0       - Π / 2
	% Joint 7       q7      0.126      0         Π / 2   */
	double UpperJointLimit[7] = { 170 * degree2Radian,120 * degree2Radian ,170 * degree2Radian ,120 * degree2Radian ,170 * degree2Radian
	,120 * degree2Radian ,175 * degree2Radian };
	double LowerJointLimit[7]= { -170 * degree2Radian,-120 * degree2Radian ,-170 * degree2Radian ,-120 * degree2Radian ,-170 * degree2Radian
		,-120 * degree2Radian ,-175 * degree2Radian };
private:
	//采样频率=FRI的通信周期
	double samplingPeriod = 0.002;

};

