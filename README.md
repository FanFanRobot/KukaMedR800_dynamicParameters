# KukaMedR800_dynamicParameters
The dynamic base parameters of KUKA Med R800, Kinematic and dynamic library. The Eigen library in C++ is necessary. 
For more details, please refer the paper "Dynamic Parameter Identification of the Kuka LBR Med Robot" with Author Fanny Ficuciello, Fan Shao. This paper will be available online after September 2023.

# Contents
We publish the base dynamic parameters of Kuka Med R800. There are two useful classes, Kinematic and DynamicModel. 

In Kinematic, there are functions about 
1) the forward kinematic of Kuka Med R800, function name "Kinematic_KukaLBRMedR800";
2) the geometry Jacobian, function name "Geometry_Jacobian_DH";
3) the PseudoInverse of the Jacobian, function name "PseudoInverseSVD_JacobianGeometry";
4) the Time Derivative of Jacobian, function name "Geometry_JacobianTimeDerivative";
5) the joint velocity, function name "JointVelocity";
6) the Cartesian velocity, function name "CartesianVelocity";
7) the orientation error between current orientation matrix and desired orientation matrix, the orientation are described by matrix, function name "OrientationError_matrix".

In DynamicModel, the base parameters of Kuka Med R800 can be found in DynamicModel.h, there are fuctions about
1. Inertia Matrix of Kuka Med R800, function name "InertiaMatrix_BaseParameters";
2. the Corriolis, Gravity and Friction terms, function name "Corriolis_Gravity_Friction".

