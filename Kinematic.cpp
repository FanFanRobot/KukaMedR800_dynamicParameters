#include "Kinematic.h"



Kinematic::Kinematic()
{
}
MatrixXd Kinematic::TransformOfAdiacentLink(double jointAngle, double alpha, double a, double d)
{
	MatrixXd T(4, 4);
	double cosJ = cos(jointAngle);
	double sinJ = sin(jointAngle);
	double cosA = cos(alpha);
	double sinA = sin(alpha);
	T << cosJ, -sinJ, 0, a,
		sinJ*cosA, cosJ*cosA, -sinA, -sinA*d,
		sinJ*sinA, cosJ*sinA, cosA, cosA*d,
		0, 0, 0, 1;

	return T;
}
MatrixXd Kinematic::Geometry_Jacobian_DH(VectorXd angles)
{
	MatrixXd J7(6, 7);
	MatrixXd T70(4, 4);
	MatrixXd T10(4, 4); MatrixXd T21(4, 4); MatrixXd T32(4, 4); MatrixXd T43(4, 4); MatrixXd T54(4, 4); MatrixXd T65(4, 4); MatrixXd T76(4, 4);
	T10 = TransformOfAdiacentLink(angles(0), alpha0, a, d1);
	T21 = TransformOfAdiacentLink(angles(1), alpha1, a, d2);
	T32 = TransformOfAdiacentLink(angles(2), alpha2, a, d3);
	T43 = TransformOfAdiacentLink(angles(3), alpha3, a, d4);
	T54 = TransformOfAdiacentLink(angles(4), alpha4, a, d5);
	T65 = TransformOfAdiacentLink(angles(5), alpha5, a, d6);
	T76 = TransformOfAdiacentLink(angles(6), alpha6, a, d7);
	T70 = T10*T21*T32*T43*T54*T65*T76;
	MatrixXd T20(4, 4); MatrixXd T30(4, 4); MatrixXd T40(4, 4); MatrixXd T50(4, 4); MatrixXd T60(4, 4);
	T60 = T10*T21*T32*T43*T54*T65;
	T50 = T10*T21*T32*T43*T54;
	T40 = T10*T21*T32*T43;
	T30 = T10*T21*T32;
	T20 = T10*T21;
	Vector3d z1(T10(0, 2), T10(1, 2), T10(2, 2)); Vector3d z2(T20(0, 2), T20(1, 2), T20(2, 2)); Vector3d z3(T30(0, 2), T30(1, 2), T30(2, 2));
	Vector3d z4(T40(0, 2), T40(1, 2), T40(2, 2)); Vector3d z5(T50(0, 2), T50(1, 2), T50(2, 2)); Vector3d z6(T60(0, 2), T60(1, 2), T60(2, 2));
	Vector3d z7(T70(0, 2), T70(1, 2), T70(2, 2));
	Vector3d p71(T70(0, 3) - T10(0, 3), T70(1, 3) - T10(1, 3), T70(2, 3) - T10(2, 3));
	Vector3d p72(T70(0, 3) - T20(0, 3), T70(1, 3) - T20(1, 3), T70(2, 3) - T20(2, 3));
	Vector3d p73(T70(0, 3) - T30(0, 3), T70(1, 3) - T30(1, 3), T70(2, 3) - T30(2, 3));
	Vector3d p74(T70(0, 3) - T40(0, 3), T70(1, 3) - T40(1, 3), T70(2, 3) - T40(2, 3));
	Vector3d p75(T70(0, 3) - T50(0, 3), T70(1, 3) - T50(1, 3), T70(2, 3) - T50(2, 3));
	Vector3d p76(T70(0, 3) - T60(0, 3), T70(1, 3) - T60(1, 3), T70(2, 3) - T60(2, 3));
	Vector3d p77(T70(0, 3) - T70(0, 3), T70(1, 3) - T70(1, 3), T70(2, 3) - T70(2, 3));

	Vector3d v1; Vector3d v2; Vector3d v3; Vector3d v4; Vector3d v5; Vector3d v6; Vector3d v7;
	v1 = z1.cross(p71); v2 = z2.cross(p72); v3 = z3.cross(p73); v4 = z4.cross(p74);
	v5 = z5.cross(p75); v6 = z6.cross(p76); v7 = z7.cross(p77);

	J7.block<3, 1>(0, 0) = v1; J7.block<3, 1>(3, 0) = z1;
	J7.block<3, 1>(0, 1) = v2; J7.block<3, 1>(3, 1) = z2;
	J7.block<3, 1>(0, 2) = v3; J7.block<3, 1>(3, 2) = z3;
	J7.block<3, 1>(0, 3) = v4; J7.block<3, 1>(3, 3) = z4;
	J7.block<3, 1>(0, 4) = v5; J7.block<3, 1>(3, 4) = z5;
	J7.block<3, 1>(0, 5) = v6; J7.block<3, 1>(3, 5) = z6;
	J7.block<3, 1>(0, 6) = v7; J7.block<3, 1>(3, 6) = z7;
	//cout << "雅克比: " << endl << J7 << endl;
	return J7;
}
MatrixXd Kinematic::PseudoInverse_JacobianGeometry(MatrixXd JacobianGeometry)
{
	MatrixXd jaInverse(7, 6);
	
	MatrixXd jaTrans;
	jaTrans = JacobianGeometry.transpose();
	jaInverse = jaTrans* ((JacobianGeometry*jaTrans).inverse());

	
	return jaInverse;
}
void Kinematic::PseudoInverseSVD_JacobianGeometry(MatrixXd Jacobian, MatrixXd* M_pinv_)
{
	bool damped = true;
	double lambda_ = damped ? 0.2 : 0.0;
	//double lambda_ = damped ? 0.05 : 0.0;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jacobian, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType sing_vals_ = svd.singularValues();
	Eigen::MatrixXd S_ = Jacobian; // copying the dimensions of M_, its content is not needed.
	S_.setZero();

	for (int i = 0; i < sing_vals_.size(); i++)
		S_(i, i) = (sing_vals_(i)) / (sing_vals_(i) * sing_vals_(i) + lambda_ * lambda_);

	*M_pinv_ = Eigen::MatrixXd(svd.matrixV() * S_.transpose() * svd.matrixU().transpose());
	//cout << *M_pinv_ << endl;
	//*M_pinv_ = (Jacobian.transpose()*Jacobian).inverse() * Jacobian.transpose();
	/**M_pinv_ = Jacobian.transpose()* (Jacobian*Jacobian.transpose()).inverse();
	cout << *M_pinv_ << endl;*/
	
}
MatrixXd Kinematic::Geometry_JacobianTimeDerivative(VectorXd joints, VectorXd jointVelocity)
{
	//t.restart();

	float q1 = joints(0); float q2 = joints(1); float q3 = joints(2); float q4 = joints(3); float q5 = joints(4); float q6 = joints(5);
	float qp1 = jointVelocity(0); float qp2 = jointVelocity(1); float qp3 = jointVelocity(2); float qp4 = jointVelocity(3); float qp5 = jointVelocity(4); float qp6 = jointVelocity(5);
	float t2 = cos(q1); float t3 = cos(q2); float t4 = cos(q3); float t5 = cos(q4); float t6 = cos(q5); float t7 = cos(q6);
	float t8 = sin(q1);
	float t9 = sin(q2);
	float t10 = sin(q3);
	float t11 = sin(q4);
	float t12 = sin(q5);
	float t13 = sin(q6);
	float t14 = d3*t3;
	float t15 = d3*t9;
	float t16 = t2*t4;
	float t17 = t3*t5;
	float t18 = t2*t10;
	float t19 = t4*t8;
	float t20 = t3*t11;
	float t21 = t5*t9;
	float t22 = t8*t10;
	float t23 = t9*t11;
	float t27 = t9*t10*t12;
	float t38 = t3*t6*t10;
	float t39 = t4*t6*t9;
	float t45 = t3*t10*t12;
	float t46 = t4*t9*t12;
	float t47 = t6*t9*t10;
	float t24 = t2*t14;
	float t25 = t8*t14;
	float t26 = t8*t23;
	float t28 = t3*t16;
	float t29 = t2*t17;
	float t30 = t4*t17;
	float t31 = t3*t18;
	float t32 = t3*t19;
	float t33 = t2*t20;
	float t34 = t2*t21;
	float t35 = t8*t17;
	float t36 = t4*t20;
	float t37 = t4*t21;
	float t40 = t3*t22;
	float t41 = t2*t23;
	float t42 = t8*t20;
	float t43 = t8*t21;
	float t44 = t4*t23;
	float t48 = d5*t10*t23;
	float t49 = t16*t21;
	float t50 = t16*t23;
	float t51 = t19*t21;
	float t52 = t6*t9*t18;
	float t53 = t6*t10*t21;
	float t55 = t19*t23;
	float t56 = t9*t12*t18;
	float t57 = t6*t9*t22;
	float t58 = t10*t12*t21;
	float t59 = t7*t10*t23;
	float t62 = t9*t12*t22;
	float t66 = -t45;
	float t67 = -t47;
	float t54 = -t28;
	float t60 = -t36;
	float t61 = -t37;
	float t63 = -t40;
	float t64 = -t41;
	float t65 = -t43;
	float t68 = -t48;
	float t69 = -t49;
	float t70 = -t51;
	float t71 = -t52;
	float t72 = -t57;
	float t73 = -t58;
	float 	t74 = -t59;
	float t75 = t18 + t32;
	float t76 = t19 + t31;
	float t77 = t17 + t44;
	float 	t78 = t23 + t30;
	float t95 = t29 + t50;
	float t98 = t35 + t55;
	float t99 = t46 + t53;
	float t79 = d5*t77;
	float t80 = t16 + t63;
	float t81 = t22 + t54;
	float t82 = t20 + t61;
	float t83 = t21 + t60;
	float 	t84 = t5*t75;
	float 	t85 = t6*t75;
	float t86 = t6*t76;
	float 	t87 = t6*t78;
	float t88 = t7*t77;
	float t89 = t11*t75;
	float t90 = t12*t75;
	float t91 = t12*t76;
	float 	t92 = t12*t78;
	float 	t93 = t13*t77;
	float 	t94 = d5*t11*t76;
	float 	t102 = t7*t11*t76;
	float t117 = t33 + t69;
	float t118 = t42 + t70;
	float t119 = t39 + t73;
	float t120 = d5*t95;
	float t121 = d5*t98;
	float t122 = t7*t95;
	float t123 = t7*t98;
	float 	t126 = t13*t99;
	float 	t96 = d5*t82;
	float 	t97 = d5*t83;
	float t100 = t5*t86;
	float t101 = t5*t91;
	float 	t103 = t6*t93;
	float 	t104 = t5*t81;
	float 	t105 = t6*t80;
	float 	t106 = t6*t81;
	float 	t107 = -t86;
	float 	t108 = t6*t82;
	float 	t109 = t7*t82;
	float 	t110 = t7*t83;
	float 	t111 = t11*t81;
	float t112 = t12*t80;
	float 	t113 = t12*t81;
	float 	t114 = t12*t82;
	float 	t115 = -t93;
	float 	t116 = d5*t11*t80;
	float 	t128 = t7*t11*t80;
	float 	t131 = t6*t117;
	float t132 = t12*t117;
	float t133 = t6*t118;
	float t135 = t12*t118;
	float 	t136 = t26 + t84;
	float 	t137 = t38 + t92;
	float 	t138 = t65 + t89;
	float 	t142 = t66 + t87;
	float t148 = -d5*(t43 - t89);
	float 	t151 = -t7*(t43 - t89);
	float 	t154 = -t13*(t45 - t87);
	float 	t155 = -t13*(t43 - t89);
	float 	t165 = t74 + t126;
	float t166 = -d7*(t59 - t126);
	float t124 = t5*t105;
	float t125 = -t100;
	float t127 = t5*t112;
	float 	t129 = -t103;
	float 	t130 = -t105;
	float t139 = t27 + t108;
	float t140 = d5*t136;
	float t141 = t34 + t111;
	float 	t143 = t6*t136;
	float t144 = t7*t136;
	float t145 = t12*t136;
	float 	t146 = t67 + t114;
	float 	t149 = t64 + t104;
	float 	t157 = -d5*(t41 - t104);
	float 	t159 = t6*t155;
	float t160 = -t6*(t41 - t104);
	float t161 = -t7*(t41 - t104);
	float t162 = -t12*(t41 - t104);
	float t163 = t56 + t131;
	float t164 = t62 + t133;
	float 	t167 = t71 + t132;
	float t168 = t72 + t135;
	float t172 = t101 + t106;
	float t176 = -d7*(t103 - t109);
	float t178 = -t13*(t100 - t113);
	float t179 = t68 + t166;
	float t184 = t110 + t154;
	float 	t134 = -t124;
	float t147 = d5*t141;
	float t150 = t7*t141;
	float 	t152 = t13*t141;
	float 	t153 = t7*t139;
	float 	t156 = t13*t139;
	float t169 = t13*t163;
	float t170 = t13*t164;
	float 	t171 = t85 + t127;
	float t174 = t113 + t125;
	float t175 = t109 + t129;
	float t180 = t112 + t143;
	float t183 = t91 + t160;
	float t185 = t130 + t145;
	float t189 = t107 + t162;
	float 	t190 = d7*t184;
	float 	t199 = t144 + t159;
	float t200 = t96 + t176;
	float 	t202 = t102 + t178;
	float t158 = t6*t152;
	float 	t173 = t90 + t134;
	float 	t181 = t88 + t156;
	float 	t186 = t115 + t153;
	float 	t187 = t7*t180;
	float t188 = t13*t180;
	float 	t191 = t7*t183;
	float 	t192 = t13*t183;
	float t195 = t122 + t169;
	float 	t196 = t123 + t170;
	float t201 = d7*t199;
	float t203 = d7*t202;
	float t210 = t97 + t190;
	float t177 = t13*t173;
	float t182 = d7*t181;
	float 	t193 = -t188;
	float t194 = -t192;
	float 	t197 = d7*t195;
	float 	t198 = d7*t196;
	float 	t204 = t158 + t161;
	float t211 = t15 + t210;
	float 	t214 = t155 + t187;
	float 	t216 = t94 + t203;
	float 	t217 = t152 + t191;
	float t218 = -d7*(t188 + t7*(t43 - t89));
	float t224 = t140 + t201;
	float t205 = t128 + t177;
	float t206 = t79 + t182;
	float 	t207 = d7*t204;
	float t212 = t120 + t197;
	float t213 = t121 + t198;
	float t215 = t151 + t193;
	float t219 = t150 + t194;
	float t226 = t148 + t218;
	float t208 = d7*t205;
	float t209 = t14 + t206;
	float t220 = t25 + t213;
	float t221 = t24 + t212;
	float t222 = d7*t219;
	float t225 = t157 + t207;
	float t223 = t116 + t208;
	float t227 = t147 + t222;
	MatrixXd dJ7(6, 7);



	dJ7.col(0) << -qp2*t220 + qp3*t223 + qp4*t224 - qp1*(t227 + t2*t15) - d7*qp6*t214 - d7*qp5*t13*(t105 - t145), qp4*(t207 - d5*(t41 - t104)) + qp3*t216 + qp2*t221 - qp1*(d7*(t188 + t7*(t43 - t89)) + t8*t15 + d5*(t43 - t89)) - d7*qp6*t217 - d7*qp5*t13*(t86 + t12*(t41 - t104)), 0.0, 0.0, 0.0, 0.0;
	dJ7.col(1) << -qp3*t2*(t48 + d7*(t59 - t126)) - qp4*t2*t200 - qp2*t2*t211 - qp1*t8*t209 - d7*qp6*t2*(t93 - t153) + d7*qp5*t2*t13*(t47 - t114), -qp3*t8*(t48 + d7*(t59 - t126)) + qp1*t2*t209 - qp4*t8*t200 - qp2*t8*t211 - d7*qp6*t8*(t93 - t153) + d7*qp5*t8*t13*(t47 - t114), -qp5*(d7*t8*t13*(t105 - t145) - d7*t2*t13*(t86 + t12*(t41 - t104))) + qp6*(d7*t2*t217 - d7*t8*t214) - qp4*(t2*(t207 - d5*(t41 - t104)) - t8*t224) - qp2*(t2*t221 + t8*t220) - qp3*(t2*t216 - t8*t223), -qp1*t2, -qp1*t8, 0.0;
	dJ7.col(2) << -qp6*(d7*t3*t214 + d7*t8*t9*(t93 - t153)) - qp1*(t3*t227 - t2*t9*t206) + qp4*(t3*t224 - t8*t9*t200) - qp5*(d7*t3*t13*(t105 - t145) - d7*t8*t9*t13*(t47 - t114)) - qp2*(t3*t213 - t9*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89)) - t3*t8*t206 + t8*t9*t210) + qp3*(t3*t223 - t8*t9*(t48 + d7*(t59 - t126))), -qp6*(d7*t3*t217 - d7*t2*t9*(t93 - t153)) + qp2*(t3*t212 - t9*t227 - t2*t3*t206 + t2*t9*t210) - qp1*(t3*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89)) - t8*t9*t206) - qp5*(d7*t3*t13*(t86 + t12*(t41 - t104)) + d7*t2*t9*t13*(t47 - t114)) + qp3*(t3*t216 + t2*t9*(t48 + d7*(t59 - t126))) + qp4*(t3*(t207 - d5*(t41 - t104)) + t2*t9*t200), qp5*(d7*t8*t9*t13*(t86 + t12*(t41 - t104)) + d7*t2*t9*t13*(t105 - t145)) + qp6*(d7*t2*t9*t214 + d7*t8*t9*t217) - qp4*(t8*t9*(t207 - d5*(t41 - t104)) + t2*t9*t224) + qp2*(t2*t3*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89)) + t2*t9*t213 - t8*t9*t212 - t3*t8*t227) - qp3*(t8*t9*t216 + t2*t9*t223), qp2*t2*t3 - qp1*t8*t9, qp1*t2*t9 + qp2*t3*t8, -qp2*t9;
	dJ7.col(3) << qp2*(t80*t210 + t3*t10*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89)) + t9*t10*t213 - t9*t22*t206) + qp6*(d7*t80*(t93 - t153) + d7*t9*t10*t214) + qp4*(t80*t200 - t9*t10*t224) + qp1*(t76*t206 + t9*t10*t227) - qp5*(d7*t13*t80*(t47 - t114) - d7*t9*t10*t13*(t105 - t145)) + qp3*(t80*(t48 + d7*(t59 - t126)) + t75*t206 + t4*t9*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89)) - t9*t10*t223), qp6*(d7*t76*(t93 - t153) + d7*t9*t10*t217) - qp1*(t80*t206 - t9*t10*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89))) + qp2*(t76*t210 - t9*t10*t212 + t9*t18*t206 - t3*t10*t227) - qp5*(d7*t13*t76*(t47 - t114) - d7*t9*t10*t13*(t86 + t12*(t41 - t104))) + qp4*(t76*t200 - t9*t10*(t207 - d5*(t41 - t104))) + qp3*(t76*(t48 + d7*(t59 - t126)) + t81*t206 - t9*t10*t216 - t4*t9*t227), qp5*(d7*t13*t76*(t105 - t145) - d7*t13*t80*(t86 + t12*(t41 - t104))) + qp6*(d7*t76*t214 - d7*t80*t217) + qp4*(t80*(t207 - d5*(t41 - t104)) - t76*t224) + qp2*(t76*t213 + t80*t212 - t9*t18*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89)) + t9*t22*t227) - qp3*(-t80*t216 + t76*t223 + t75*t227 + t81*(d7*(t188 + t7*(t43 - t89)) + d5*(t43 - t89))), qp1*t80 - qp3*t81 - qp2*t9*t18, qp1*t76 + qp3*t75 - qp2*t9*t22, -qp2*t3*t10 - qp3*t4*t9;
	dJ7.col(4) << qp3*(t77*t208 - d7*(t43 - t89)*(t59 - t126) - t11*t80*t182 + d7*t10*t23*(t188 + t7*(t43 - t89))) + qp4*(t77*t201 - t136*t182 + d7*(t43 - t89)*(t103 - t109) + d7*t82*(t188 + t7*(t43 - t89))) - qp2*(t77*t198 - t98*t182 + t83*t218 + t190*(t43 - t89)) - qp1*(t77*t222 - t141*t182) - qp5*(d7*t155*(t47 - t114) + d7*t93*(t105 - t145)) - qp6*(d7*(t43 - t89)*(t93 - t153) + d7*t77*t214), qp4*(t77*t207 - t82*t222 + t182*(t41 - t104) - d7*t141*(t103 - t109)) + qp2*(t77*t197 - t95*t182 - t83*t222 + t141*t190) - qp5*(d7*t93*(t86 + t12*(t41 - t104)) + d7*t152*(t47 - t114)) - qp6*(d7*t77*t217 - d7*t141*(t93 - t153)) + qp3*(t77*t203 - t10*t23*t222 - t11*t76*t182 + d7*t141*(t59 - t126)) + qp1*(t77*t218 + t182*(t43 - t89)), qp5*(d7*t152*(t105 - t145) + d7*t13*(t86 + t12*(t41 - t104))*(t43 - t89)) - qp3*(t141*t208 + t203*(t43 - t89) + t11*t76*t218 - t11*t80*t222) + qp6*(d7*t141*t214 + d7*t217*(t43 - t89)) - qp2*(t98*t222 - t141*t198 + t197*(t43 - t89) - d7*t95*(t188 + t7*(t43 - t89))) - qp4*(t141*t201 - t136*t222 + t207*(t43 - t89) + d7*(t188 + t7*(t43 - t89))*(t41 - t104)), qp2*t95 - qp1*(t43 - t89) - qp4*(t41 - t104) + qp3*t11*t76, qp2*t98 - qp4*t136 + qp1*t141 - qp3*t11*t80, -qp2*t83 - qp4*t82 - qp3*t10*t23;
	dJ7.col(5) << -qp6*(d7*(t93 - t153)*(t105 - t145) + d7*t214*(t47 - t114)) - qp2*(t198*(t47 - t114) - t182*(t57 - t135) + t190*(t105 - t145) + d7*t137*(t188 + t7*(t43 - t89))) - qp3*(t171*t182 - t208*(t47 - t114) + d7*(t59 - t126)*(t105 - t145) + d7*t119*(t188 + t7*(t43 - t89))) + qp4*(t201*(t47 - t114) + d7*(t103 - t109)*(t105 - t145) - t12*t182*(t43 - t89) + d7*t12*t77*(t188 + t7*(t43 - t89))) - qp5*(t139*t218 + t180*t182) - qp1*(t182*(t86 + t12*(t41 - t104)) + t222*(t47 - t114)), qp1*(t218*(t47 - t114) + t182*(t105 - t145)) + qp3*(t166*(t86 + t12*(t41 - t104)) + t119*t222 - t172*t182 + t203*(t47 - t114)) - qp2*(t190*(t86 + t12*(t41 - t104)) - t137*t222 - t197*(t47 - t114) + t182*(t52 - t132)) - qp5*(t139*t222 + t182*t183) - qp6*(d7*(t86 + t12*(t41 - t104))*(t93 - t153) + d7*t217*(t47 - t114)) + qp4*(t207*(t47 - t114) + d7*(t86 + t12*(t41 - t104))*(t103 - t109) - t12*t77*t222 + t12*t141*t182), qp4*(t201*(t86 + t12*(t41 - t104)) - t207*(t105 - t145) + t12*t141*t218 + t12*t222*(t43 - t89)) - qp6*(d7*t214*(t86 + t12*(t41 - t104)) - d7*t217*(t105 - t145)) + qp3*(t208*(t86 + t12*(t41 - t104)) + t171*t222 - t203*(t105 - t145) + d7*t172*(t188 + t7*(t43 - t89))) - qp2*(t198*(t86 + t12*(t41 - t104)) + t218*(t52 - t132) + t222*(t57 - t135) + t197*(t105 - t145)) + qp5*(t180*t222 + d7*t183*(t188 + t7*(t43 - t89))), qp3*t172 + qp5*t183 + qp2*(t52 - t132) - qp1*(t105 - t145) - qp4*t12*t141, -qp1*(t86 + t12*(t41 - t104)) - qp3*t171 - qp5*t180 + qp2*(t57 - t135) - qp4*t12*(t43 - t89), qp3*t119 + qp2*t137 - qp5*t139 - qp4*t12*t77;
	dJ7.col(6) << 0.0, 0.0, 0.0, -qp1*(t188 + t7*(t43 - t89)) + qp2*t195 + qp3*t202 + qp4*t204 - qp6*t217 - qp5*t13*(t86 + t12*(t41 - t104)), qp2*t196 - qp4*t199 - qp3*t205 + qp1*t219 + qp6*t214 + qp5*t13*(t105 - t145), -qp2*t184 - qp3*(t59 - t126) + qp4*(t103 - t109) - qp6*(t93 - t153) + qp5*t13*(t47 - t114);

	//cout <<"dJ7: "<< dJ7 << endl;

	/*double ela = t.elapsed() * 1000;
	cout << ela << endl;*/
	return dJ7;
}
MatrixXd Kinematic::Kinematic_KukaLBRMedR800(double angles[])
{
	MatrixXd T70(4, 4);
	MatrixXd T10(4, 4); MatrixXd T21(4, 4); MatrixXd T32(4, 4); MatrixXd T43(4, 4); MatrixXd T54(4, 4); MatrixXd T65(4, 4); MatrixXd T76(4, 4);
	T10 = TransformOfAdiacentLink(angles[0], alpha0, a, d1);
	T21 = TransformOfAdiacentLink(angles[1], alpha1, a, d2);
	T32 = TransformOfAdiacentLink(angles[2], alpha2, a, d3);
	T43 = TransformOfAdiacentLink(angles[3], alpha3, a, d4);
	T54 = TransformOfAdiacentLink(angles[4], alpha4, a, d5);
	T65 = TransformOfAdiacentLink(angles[5], alpha5, a, d6);
	T76 = TransformOfAdiacentLink(angles[6], alpha6, a, d7);
	T70 = T10*T21*T32*T43*T54*T65*T76;

	return T70;
}
MatrixXd Kinematic::Kinematic_KukaLBRMedR800(VectorXd angles)
{
	MatrixXd T70(4, 4);
	MatrixXd T10(4, 4); MatrixXd T21(4, 4); MatrixXd T32(4, 4); MatrixXd T43(4, 4); MatrixXd T54(4, 4); MatrixXd T65(4, 4); MatrixXd T76(4, 4);
	T10 = TransformOfAdiacentLink(angles(0), alpha0, a, d1);
	T21 = TransformOfAdiacentLink(angles(1), alpha1, a, d2);
	T32 = TransformOfAdiacentLink(angles(2), alpha2, a, d3);
	T43 = TransformOfAdiacentLink(angles(3), alpha3, a, d4);
	T54 = TransformOfAdiacentLink(angles(4), alpha4, a, d5);
	T65 = TransformOfAdiacentLink(angles(5), alpha5, a, d6);
	T76 = TransformOfAdiacentLink(angles(6), alpha6, a, d7);
	T70 = T10*T21*T32*T43*T54*T65*T76;

	return T70;
}
/*
Param 1: dj,
it is the difference of the joint position for the current time period,
Param 2: dj_,
it is the difference of the joint position for the last time peroid.
Velocity calculation from the method "2nd Fixed-Time TSE(TAYLOR SERIES EXPANSION) Estimators",
please refer paper "Analysis of Algorithms for Velocity Estimation from Discrete Position Versus Time Data" for more details.
*/
VectorXd Kinematic::JointVelocity(VectorXd dj, VectorXd dj_)
{
	VectorXd vj_now(6);
	vj_now = dj / samplingPeriod;
	/*VectorXd vj_last(6);
	vj_last = dj_ / samplingPeriod;
	VectorXd aj_now(6);
	aj_now = (vj_now - vj_last) / samplingPeriod;

	VectorXd vj(6);
	vj = vj_now + 0.5f*aj_now;*/


	return vj_now;
}

/*
Calculate the CartesianVelocity,
Param 1: vj,
it is joint velocity (rad/s),
Param 2: Jacobian_geometry,
it is the geometry jacobian (6*7)
*/
VectorXd Kinematic::CartesianVelocity(VectorXd vj, MatrixXd Jacobian_geometry)
{
	VectorXd cartV(6);
	cartV = Jacobian_geometry*vj;
	return cartV;
}
/*
Calculate the orientation error between current orientation matrix and desired orientation matrix. The orientation are described by matrix.
Param 1: desiredPose,
it is the desired matrix expressed in the base frame. The current is also expressed in the base frame.
Param 2: currentPose,
the current matrix, also expressed in the base frame.
More details please refer to paper "Pose Control of Robot Manipulators Using Different Orientation Representations: A Comparative Review".
*/
VectorXd Kinematic::OrientationError_matrix(MatrixXd desiredPose, MatrixXd currentPose)
{
	MatrixXd oriError(4, 4);
	oriError = desiredPose*(currentPose.transpose());
	Vector3d oriErrVec(0.5f*(oriError(2, 1) - oriError(1, 2)), 0.5f*(oriError(0, 2) - oriError(2, 0)), 0.5f*(oriError(1, 0) - oriError(0, 1)));
	return oriErrVec;
}
/*
这个旋转表示的方式可能还有错误或者问题，待定
*/
Vector3d Kinematic::OrientationError_quaternion(const Quaterniond &orientation_d, Quaterniond orientation)
{
	// Orientation error
	if (orientation_d.coeffs().dot(orientation.coeffs()) < 0.0)
	{
		orientation.coeffs() << -orientation.coeffs();
	}
	// "difference" quaternion
	const Quaterniond error_quaternion(orientation * orientation_d.inverse());
	// convert to axis angle
	AngleAxisd error_quaternion_angle_axis(error_quaternion);
	return error_quaternion_angle_axis.axis() * error_quaternion_angle_axis.angle();
}
/*
Cartesian Position error. (6*1 vector)
Param 1:desiredPose,
it is the desired matrix expressed in the base frame. The current is also expressed in the base frame.
Param 2: currentPose,
the current matrix, also expressed in the base frame.
*/
VectorXd Kinematic::PosiError(MatrixXd desiredPose, MatrixXd currentPose)
{
	Vector3d oriErrVec = OrientationError_matrix(desiredPose, currentPose);
	Matrix3d oriC;
	oriC << currentPose(0, 0), currentPose(0, 1), currentPose(0, 2),
		currentPose(1, 0), currentPose(1, 1), currentPose(1, 2),
		currentPose(2, 0), currentPose(2, 1), currentPose(2, 2);
	/*Quaterniond orientation_c{ Eigen::Quaterniond::Identity() };
	orientation_c = oriC;
	Vector3d oriErrVec = OrientationError_quaternion(this->orientation_d_, orientation_c);*/
	VectorXd posiErr(6);
	posiErr(0) = desiredPose(0, 3) - currentPose(0, 3);
	posiErr(1) = desiredPose(1, 3) - currentPose(1, 3);
	posiErr(2) = desiredPose(2, 3) - currentPose(2, 3);
	posiErr(3) = oriErrVec(0);
	posiErr(4) = oriErrVec(1);
	posiErr(5) = oriErrVec(2);

	//test-s
	/*cout << desiredPose << endl;
	cout << currentPose << endl;
	cout << posiErr << endl;*/
	//test-e
	return posiErr;
}

VectorXd Kinematic::VelError(VectorXd desiredVel, VectorXd currentVel)
{
	VectorXd velErr(6);
	velErr = desiredVel - currentVel;
	return velErr;
}
Kinematic::~Kinematic()
{
}
