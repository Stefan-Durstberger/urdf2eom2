//
//  SymbolicForwardDynamics.hpp
//  urdf2eom
//
//  Created by Stefan on 06.08.18.
//  Copyright © 2018 Stefan Durstberger. All rights reserved.
//

#ifndef SymbolicForwardDynamics_hpp
#define SymbolicForwardDynamics_hpp

#include <stdio.h>

#include <iostream>
#include <limits>
#include <assert.h>
#include <string.h>

#include "srbdl/srbdl_mathutils.h"
#include "srbdl/Logging.h"

#include "srbdl/Model.h"
#include "srbdl/Joint.h"
#include "srbdl/Body.h"
#include "srbdl/Dynamics.h"
#include "srbdl/Kinematics.h"

#include "symbolicc++/symbolicc++.h"

#include "SymbolicMath.hpp"

using namespace SymbolicRigidBodyDynamics;
using namespace SymbolicMath;
using namespace std;


static inline void jcalc ( Model &model, unsigned int joint_id );
static inline void SymbolicForwardDynamics ( Model &model, std::vector<SymSpatialVector> *f_ext = NULL );


// State information
/// \brief The spatial velocity of the bodies
std::vector<SymSpatialVector> sym_v;
/// \brief The spatial acceleration of the bodies
std::vector<SymSpatialVector> sym_a;

////////////////////////////////////
// Joints

/// \brief All joints

//std::vector<Joint> mJoints;
/// \brief The joint axis for joint i
std::vector<SymSpatialVector> sym_S;

// Joint state variables
std::vector<SymSpatialTransform> sym_X_J;
std::vector<SymSpatialVector> sym_v_J;
std::vector<SymSpatialVector> sym_c_J;

//std::vector<unsigned int> mJointUpdateOrder;

/// \brief Transformations from the parent body to the frame of the joint.
// It is expressed in the coordinate frame of the parent.
//std::vector<Math::SpatialTransform> X_T;
/// \brief The number of fixed joints that have been declared before
///  each joint.
//std::vector<unsigned int> mFixedJointCount;

////////////////////////////////////
// Special variables for joints with 3 degrees of freedom
/// \brief Motion subspace for joints with 3 degrees of freedom
//std::vector<Math::Matrix63> multdof3_S;
//std::vector<Math::Matrix63> multdof3_U;
//std::vector<Math::Matrix3d> multdof3_Dinv;
//std::vector<Math::Vector3d> multdof3_u;
//std::vector<unsigned int> multdof3_w_index;
//
//std::vector<CustomJoint*> mCustomJoints;

////////////////////////////////////
// Dynamics variables

/// \brief The velocity dependent spatial acceleration
std::vector<SymSpatialVector> sym_c;
/// \brief The spatial inertia of the bodies
std::vector<SymSpatialMatrix> sym_IA;
/// \brief The spatial bias force
std::vector<SymSpatialVector> sym_pA;
/// \brief Temporary variable U_i (RBDA p. 130)
std::vector<SymSpatialVector> sym_U;
/// \brief Temporary variable D_i (RBDA p. 130)
Symbolic d;
/// \brief Temporary variable u (RBDA p. 130)
Symbolic u;
/// \brief Internal forces on the body (used only InverseDynamics())
//std::vector<Math::SpatialVector> f;
/// \brief The spatial inertia of body i (used only in
///  CompositeRigidBodyAlgorithm())
std::vector<SymSpatialRigidBodyInertia> sym_I;
//std::vector<SymSpatialRigidBodyInertia> Ic;
//std::vector<Math::SpatialVector> hc;
//std::vector<Math::SpatialVector> hdotc;

////////////////////////////////////
// Bodies

/** \brief Transformation from the parent body to the current body
 * \f[
 *	X_{\lambda(i)} = {}^{i} X_{\lambda(i)}
 * \f]
 */
std::vector<SymSpatialTransform> sym_X_lambda;
/// \brief Transformation from the base to bodies reference frame
std::vector<SymSpatialTransform> sym_X_base;
std::vector<SymSpatialMatrix> IA;

/// \brief All bodies that are attached to a body via a fixed joint.
//std::vector<FixedBody> mFixedBodies;


static inline void SymbolicForwardDynamics ( Model &model, std::vector<SymSpatialVector> *f_ext ) {
	
	for (unsigned int i = 0; i < model.v_J.size(); i++) {
		list<Symbolic> symList;
		for (unsigned int j = 0; j < model.v_J[i].size(); j++) {
			symList.emplace_back(Symbolic (model.v_J[i][j] ));
		}
		sym_v_J.push_back(SymSpatialVector (symList));
	}
	for (unsigned int i = 0; i < model.S.size(); i++) {
		list<Symbolic> symList;
		for (unsigned int j = 0; j < model.S[i].size(); j++) {
			symList.emplace_back(Symbolic (model.S[i][j] ));
		}
		sym_S.push_back(SymSpatialVector (symList));
	}
	for (unsigned int i = 0; i < model.I.size(); i++) {
		Symbolic sym_copy_h (list<Symbolic> {Symbolic (model.I[i].h[0]), Symbolic(model.I[i].h[1]), Symbolic(model.I[i].h[2]) });
		Symbolic sym_copy_I (list<list<Symbolic>> {
			{Symbolic (model.I[i].Ixx), symZero, 				  symZero},
			{Symbolic (model.I[i].Iyx), Symbolic(model.I[i].Iyy), symZero},
			{Symbolic (model.I[i].Izx), Symbolic(model.I[i].Izy), Symbolic(model.I[i].Izz)}
			});
		SymSpatialRigidBodyInertia sym_copy_rbi(Symbolic (model.I[i].m), sym_copy_h, sym_copy_I );
		sym_I.emplace_back(sym_copy_rbi);
	}
	
//	clog << model.X_J.size() << endl;
//	clog << model.X_lambda.size() << endl;
//	clog << model.X_base.size()  << endl;
//	clog << model.v.size()  << endl;
//	clog << model.c_J.size()  << endl;
//	clog << model.c.size()  << endl;
//	clog << model.IA.size()  << endl;
//	clog << model.pA.size()  << endl;
//	clog << model.U.size()  << endl;
//	cout << "dof = " << model.dof_count << endl;
	
	for (unsigned int i = 0; i < model.dof_count + 1; i++) {
		sym_X_J.push_back 		(SymSpatialTransform());
		sym_X_lambda.push_back	(SymSpatialTransform());
		sym_X_base.push_back 	(SymSpatialTransform());
		sym_v.push_back 		(symSpatialVectorZero);
		sym_c_J.push_back 		(symSpatialVectorZero);
		sym_c.push_back 		(symSpatialVectorZero);
		sym_IA.push_back		(symSpatialMatrixIdentity);
		sym_pA.push_back		(symSpatialVectorZero);
		sym_U.push_back			(symSpatialVectorZero);
		//sym_Ic.push_back		(symSRBIZero);
	}
	
	//clog << symSpatialVectorZero.v << endl;
	
	Math::SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);
	
	unsigned int i = 0;
	
	// Reset the velocity of the root body
	model.v[0].setZero();
	
//	//Test SymSpatialTransform.apply
//	clog << endl;
//	clog << "---------------------------------------" << endl;
//	clog << "[Test START]  SymSpatialTransform.apply" << endl;
//	clog << endl;
//	SymSpatialVector testSpatVec;
//	for (unsigned int i = 0; i < 6; i++) {
//		testSpatVec.v(i) = Symbolic ((int)i);
//	}
//	clog << "testSpatVec = " << testSpatVec.v << endl;
//
//	SymSpatialTransform testSpatTransform;
//	for (unsigned int i = 0; i < 3; i++) {
//		for (unsigned int j = 0; j < 3; j++) {
//			testSpatTransform.E(i,j) = Symbolic ( (int)(i*j) );
//		}
//		testSpatTransform.r(i) = Symbolic ((int)i);
//	}
//	clog << "testSpatTransform.E = " << testSpatTransform.E << endl;
//	clog << "testSpatTransform.r = " << testSpatTransform.r << endl;
//
//	SymSpatialVector testSpatVec2;
//	testSpatVec2 = testSpatTransform.apply(testSpatVec);
//	clog << "testSpatVec2.v = " << testSpatVec2.v << endl;
//	clog << "[Test END]  SymSpatialTransform.apply" << endl;
//	clog << "-------------------------------------" << endl;
//	clog << endl;
//	// Test SymSpatialTransform.apply
	
	for (i = 1; i < model.mBodies.size(); i++) {
		// START_TIMING 1
		auto start1 = std::chrono::high_resolution_clock::now();
		// START_TIMING 1
		
		clog << "i = " << i << endl;
		
		unsigned int lambda = model.lambda[i];
		
		jcalc (model, i);
		
		if (lambda != 0) {
			sym_X_base[i] = sym_X_lambda[i] * SpatialTrans2SymSpatialTrans (model.X_base[lambda]);
		}
		else
			sym_X_base[i] = sym_X_lambda[i];

//		clog << "sym_X_base[" << i << "].E = " << sym_X_base[i].E << endl;
//		clog << "sym_X_base[" << i << "].r = " << sym_X_base[i].r << endl;
//		clog << "sym_X_lambda[" << i << "].E = " << sym_X_lambda[i].E << endl;
//		clog << "sym_X_lambda[" << i << "].r = " << sym_X_lambda[i].r << endl;
//		clog << "sym_v[" << lambda << "].v = " << sym_v[lambda].v << endl;
//		clog << "sym_v_J[" << i << "].v = " << sym_v_J[i].v << endl;
		
		sym_v[i] = sym_X_lambda[i].apply( sym_v[lambda]) + sym_v_J[i];
		
//		clog << "sym_v[" << i << "].v = " << sym_v[i].v << endl;
		
		sym_c[i] = sym_c_J[i] + crossm(sym_v[i],sym_v_J[i]);
		
//		clog << "sym_c[" << i << "].v = " << sym_c[i].v << endl;
		
		sym_I[i].setSpatialMatrix (sym_IA[i].M);
		
		//clog << "sym_c[" << i << "].v = " << sym_c[i].v << endl;

		sym_pA[i] = crossf(sym_v[i].v,sym_I[i] * sym_v[i]);
		
		if (f_ext != NULL && (*f_ext)[i] != symSpatialVectorZero) {
			clog << "External force (" << i << ") = " << sym_X_base[i].toMatrixAdjoint() * (*f_ext)[i] << std::endl;
			sym_pA[i] -= sym_X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
		
		//std::clog << "sym_pA[" << i << "].v = " << sym_pA[i].v << std::endl;
		
		// END_TIMING 1
		auto finish1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = finish1 - start1;
		std::clog << "[TIME 1] Time elapsed = "<< elapsed1.count() << std::endl;
		// END_TIMING 1
	}
	
//	for (i = static_cast<unsigned int>(model.mBodies.size()) - 1; i > 0; i--) {
//		unsigned int q_index = model.mJoints[i].q_index;
//
//		if (model.mJoints[i].mDoFCount == 1
//			&& model.mJoints[i].mJointType != JointTypeCustom) {
//
//			sym_U[i] = sym_IA[i] * sym_S[i];
//			model.d[i] = model.S[i].dot(model.U[i]);
//			model.u[i] = Tau[q_index] - model.S[i].dot(model.pA[i]);
//			//      LOG << "u[" << i << "] = " << model.u[i] << std::endl;
//
//			unsigned int lambda = model.lambda[i];
//			if (lambda != 0) {
//				SpatialMatrix Ia =    model.IA[i]
//				- model.U[i]
//				* (model.U[i] / model.d[i]).transpose();
//
//				SpatialVector pa =  model.pA[i]
//				+ Ia * model.c[i]
//				+ model.U[i] * model.u[i] / model.d[i];
//
//#ifdef EIGEN_CORE_H
//				model.IA[lambda].noalias()
//				+= model.X_lambda[i].toMatrixTranspose()
//				* Ia * model.X_lambda[i].toMatrix();
//				model.pA[lambda].noalias()
//				+= model.X_lambda[i].applyTranspose(pa);
//#else
//				model.IA[lambda]
//				+= model.X_lambda[i].toMatrixTranspose()
//				* Ia * model.X_lambda[i].toMatrix();
//
//				model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
//#endif
//				LOG << "pA[" << lambda << "] = "
//				<< model.pA[lambda].transpose() << std::endl;
//			}
//		} else if (model.mJoints[i].mDoFCount == 3
//				   && model.mJoints[i].mJointType != JointTypeCustom) {
//			model.multdof3_U[i] = model.IA[i] * model.multdof3_S[i];
//#ifdef EIGEN_CORE_H
//			model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose()
//									  * model.multdof3_U[i]).inverse().eval();
//#else
//			model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose()
//									  * model.multdof3_U[i]).inverse();
//#endif
//			Vector3d tau_temp(Tau.block(q_index,0,3,1));
//			model.multdof3_u[i] = tau_temp
//			- model.multdof3_S[i].transpose() * model.pA[i];
//
//			// LOG << "multdof3_u[" << i << "] = "
//			//                      << model.multdof3_u[i].transpose() << std::endl;
//			unsigned int lambda = model.lambda[i];
//			if (lambda != 0) {
//				SpatialMatrix Ia = model.IA[i]
//				- model.multdof3_U[i]
//				* model.multdof3_Dinv[i]
//				* model.multdof3_U[i].transpose();
//				SpatialVector pa = model.pA[i]
//				+ Ia
//				* model.c[i]
//				+ model.multdof3_U[i]
//				* model.multdof3_Dinv[i]
//				* model.multdof3_u[i];
//#ifdef EIGEN_CORE_H
//				model.IA[lambda].noalias()
//				+= model.X_lambda[i].toMatrixTranspose()
//				* Ia
//				* model.X_lambda[i].toMatrix();
//
//				model.pA[lambda].noalias()
//				+= model.X_lambda[i].applyTranspose(pa);
//#else
//				model.IA[lambda]
//				+= model.X_lambda[i].toMatrixTranspose()
//				* Ia
//				* model.X_lambda[i].toMatrix();
//
//				model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
//#endif
//				LOG << "pA[" << lambda << "] = "
//				<< model.pA[lambda].transpose()
//				<< std::endl;
//			}
//		} else if (model.mJoints[i].mJointType == JointTypeCustom) {
//			unsigned int kI   = model.mJoints[i].custom_joint_index;
//			unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;
//			model.mCustomJoints[kI]->U =
//			model.IA[i] * model.mCustomJoints[kI]->S;
//
//#ifdef EIGEN_CORE_H
//			model.mCustomJoints[kI]->Dinv
//			= (model.mCustomJoints[kI]->S.transpose()
//			   * model.mCustomJoints[kI]->U).inverse().eval();
//#else
//			model.mCustomJoints[kI]->Dinv
//			= (model.mCustomJoints[kI]->S.transpose()
//			   * model.mCustomJoints[kI]->U).inverse();
//#endif
//			VectorNd tau_temp(Tau.block(q_index,0,dofI,1));
//			model.mCustomJoints[kI]->u = tau_temp
//			- model.mCustomJoints[kI]->S.transpose() * model.pA[i];
//
//			//      LOG << "multdof3_u[" << i << "] = "
//			//      << model.multdof3_u[i].transpose() << std::endl;
//			unsigned int lambda = model.lambda[i];
//			if (lambda != 0) {
//				SpatialMatrix Ia = model.IA[i]
//				- (model.mCustomJoints[kI]->U
//				   * model.mCustomJoints[kI]->Dinv
//				   * model.mCustomJoints[kI]->U.transpose());
//				SpatialVector pa =  model.pA[i]
//				+ Ia * model.c[i]
//				+ (model.mCustomJoints[kI]->U
//				   * model.mCustomJoints[kI]->Dinv
//				   * model.mCustomJoints[kI]->u);
//
//#ifdef EIGEN_CORE_H
//				model.IA[lambda].noalias() += model.X_lambda[i].toMatrixTranspose()
//				* Ia
//				* model.X_lambda[i].toMatrix();
//				model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
//#else
//				model.IA[lambda] += model.X_lambda[i].toMatrixTranspose()
//				* Ia
//				* model.X_lambda[i].toMatrix();
//				model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
//#endif
//				LOG << "pA[" << lambda << "] = "
//				<< model.pA[lambda].transpose()
//				<< std::endl;
//			}
//		}
//	}
}


static inline void jcalc ( Model &model, unsigned int joint_id ) {
	
	// exception if we calculate it for the root body
	assert (joint_id > 0);
	
	Symbolic q("q", model.dof_count);
	Symbolic qdot("qdot", model.dof_count);
	
	if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
		//		std::clog << "model.v_J.size() = " << model.v_J.size() << std::endl;
		//		std::clog << "model.X_J.size() = " << model.X_J.size() << std::endl;
		//		std::clog << "sym_X_J.size() = " << sym_X_J.size() << std::endl;
		//		std::clog << "joint_id = " << joint_id << std::endl;
		//		std::clog << "sym_v_J = " << sym_v_J << std::endl;
		
		sym_X_J[joint_id] = SymXrotx (q(model.mJoints[joint_id].q_index));
		sym_v_J[joint_id].v(0) = qdot(model.mJoints[joint_id].q_index);
		//		std::clog << "[jcalc]: sym_X_J = " << std::endl;
		//		std::clog << sym_X_J[joint_id].E << std::endl;
		//		std::clog << sym_X_J[joint_id].r << std::endl;
	} else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
		sym_X_J[joint_id] = SymXroty (q(model.mJoints[joint_id].q_index));
		sym_v_J[joint_id].v(1) = qdot(model.mJoints[joint_id].q_index);
	} else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
		sym_X_J[joint_id] = SymXrotz (q(model.mJoints[joint_id].q_index));
		sym_v_J[joint_id].v(2) = qdot(model.mJoints[joint_id].q_index);
		/* Joint Helical, Sherical, EulerZYX, EulerXYZ, EulerYXZ, TranlationalXYZ, Custom not available */
	} else {
		std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
		abort();
	}
	sym_X_lambda[joint_id] = sym_X_J[joint_id] * SpatialTrans2SymSpatialTrans(model.X_T[joint_id]);
}


#endif /* SymbolicForwardDynamics_hpp */
