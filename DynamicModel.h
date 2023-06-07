#pragma once
#include <math.h>
#include <iostream>
#include "mkl.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include<boost/progress.hpp>
#include<boost/timer.hpp>
#include "Kinematic.h"
#define _USE_MATH_DEFINES
using namespace Eigen;
using namespace std;
class DynamicModel
{
public:
	DynamicModel(Kinematic *k);
	~DynamicModel();
	MatrixXd InertiaMatrix_BaseParameters(VectorXd joints);
	MatrixXd InertiaWeightedPseudoInverse_JacobianGeometry(MatrixXd JacobianGeometry, MatrixXd Inertia);
	VectorXd Corriolis_Gravity_Friction(VectorXd joints, VectorXd jointsVelocity);
	//base parameters

	float Beta1 = 0.0035;	float Beta2 = 0.0007;	float Beta3 = 0.1094;	float Beta4 = -0.0433;	float Beta5 = -0.0252;
	float Beta6 = 0.0093;	float Beta7 = -0.0269;	float Beta8 = -0.0542;	float Beta9 = 0.0109;	float Beta10 = 2.2794;
	float Beta11 = 0.0779;	float Beta12 = 0.0453;	float Beta13 = 0.1772; float Beta14 = 0.0579;	float Beta15 = -0.0402;
	float Beta16 = 0.1114;	float Beta17 = 0.0093;	float Beta18 = -0.1248;	float Beta19 = 5.5162;
	float Beta20 = 0.1181;	float Beta21 = 0.0562;	float Beta22 = -0.1613;	float Beta23 = -0.1806;
	float Beta24 = -0.1801;	float Beta25 = 0.1593;	float Beta26 = -0.1112;	float Beta27 = -0.0708;
	float Beta28 = -0.0701;	float Beta29 = -0.0399;	float Beta30 = -0.5550;	float Beta31 = 0.0488;
	float Beta32 = -0.5981;	float Beta33 = 0.1813;	float Beta34 = -3.4424;	float Beta35 = 0.1860;
	float Beta36 = 0.1135;	float Beta37 = -3.5834;	float Beta38 = 0.8261;	float Beta39 = 0.5910;
	float Beta40 = -7.6326;	float Beta41 = -10.8094;	float Beta49 = 0.0190;	float Beta50 = 0.2818;
	float Beta51 = 0.1106;	float Beta52 = 0.0165;	float Beta53 = 1.7202;	float Beta54 = 0.0221;
	float Beta55 = 0.1821;
	float Beta42 = -0.2507; float Beta43 = -0.1182; float Beta44 = -0.1484; float Beta45 = -0.1445; float Beta46 = -0.0048;
	float Beta47 = -0.4559; float Beta48 = 0.0584;
	float Beta56 = -0.2481;
	float Beta57 = -1.4887;
	float Beta58 = -0.1271;
	float Beta59 = -0.2397;
	float Beta60 = -0.2916;
	float Beta61 = -0.9850;
	float Beta62 = -0.0613;
private:
	double alpha0, alpha1, alpha2, alpha3,
		alpha4, alpha5, alpha6;
	double a;
	double d1, d2, d3, d4, d5, d6, d7;
	Kinematic *kine;
};
