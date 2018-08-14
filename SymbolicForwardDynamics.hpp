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
#include "SymbolicDataTypes.hpp"
#include "HelperFunctions.hpp"

using namespace SymbolicRigidBodyDynamics;
using namespace SymbolicMath;
using namespace SymbolicDataTypes;
using namespace std;


static inline void jcalc ( Model &model, unsigned int joint_id );
static inline void SymbolicForwardDynamics ( Model &model, std::vector<SymSpatialVector> *f_ext = NULL );

Timer timer1;

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
std::vector<Symbolic> sym_D;
/// \brief Temporary variable u (RBDA p. 130)
std::vector<Symbolic> sym_u;
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
		sym_a.push_back			(symSpatialVectorZero);
	}
	
	for (unsigned int i = 0; i < model.mBodies.size(); i++) {
		sym_u.push_back(symZero);
		sym_D.push_back(symZero);
	}
	
	Symbolic Tau("tau", model.dof_count);
	
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
	
	SymSpatialVector spatial_gravity (list<Symbolic> {symZero, symZero, symZero, Symbolic (model.gravity[0]), Symbolic (model.gravity[1]), Symbolic (model.gravity[2])} );
	
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
	
	/* FIRST LOOP */
	for (i = 1; i < model.mBodies.size(); i++) {
		clog << "[LOOP 1] Iteration "<< i << "/" << model.mBodies.size() - 1 << endl << endl;
		
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
//		clog << "sym_I[" << i << "].m = " << sym_I[i].m << endl;
//		clog << "sym_I[" << i << "].h = " << sym_I[i].h << endl;
//		clog << "sym_I[" << i << "].Ixx = " << sym_I[i].Ixx << endl;
//		clog << "sym_I[" << i << "].Iyx = " << sym_I[i].Iyx << ";    Iyy = " << sym_I[i].Iyy << endl;
//		clog << "sym_I[" << i << "].Izx = " << sym_I[i].Izx << ";    Izy = " << sym_I[i].Izy << ";    Izz = " << sym_I[i].Izz << endl;
		
		sym_v[i] = sym_X_lambda[i].apply (sym_v[lambda]) + sym_v_J[i];
		
//		clog << "sym_v[" << i << "].v = " << sym_v[i].v << endl;
		
		sym_c[i] = sym_c_J[i] + crossm(sym_v[i],sym_v_J[i]);
		
//		clog << "sym_c[" << i << "].v = " << sym_c[i].v << endl;
		
		sym_I[i].setSpatialMatrix (sym_IA[i]);

//		clog << "sym_IA[" << i << "] = " << sym_IA[i].M << std::endl;

		sym_pA[i] = crossf (sym_v[i],sym_I[i] * sym_v[i]);
		
		if (f_ext != NULL && (*f_ext)[i] != symSpatialVectorZero) {
			clog << "External force (" << i << ") = " << sym_X_base[i].toMatrixAdjoint() * (*f_ext)[i] << std::endl;
			sym_pA[i] -= sym_X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
		
//		std::clog << "sym_pA[" << i << "].v = " << sym_pA[i].v << std::endl;
	} /* FIRST LOOP */
	
	//Evaluate Model:
	//cout << "sym_pA[" << i - 1 << "] = " << sym_pA[i-1].v << endl;
	
	/* SECOND LOOP */
	for (i = static_cast<unsigned int>(model.mBodies.size()) - 1; i > 0; i--) {
		clog << "[LOOP 2] Iteration "<< model.mBodies.size() - i << "/" << model.mBodies.size() - 1 << endl << endl;
		timer1.start();
		
		unsigned int q_index = model.mJoints[i].q_index;
		
		if (model.mJoints[i].mDoFCount == 1
			&& model.mJoints[i].mJointType != JointTypeCustom) {
			
//			clog << "sym_S[" << i << "] = " << sym_S[i].v << std::endl;

			sym_U[i] = sym_IA[i] * sym_S[i];
			
//			clog << "sym_U[" << i << "] = " << sym_U[i].v << std::endl;
			
			sym_D[i] = sym_S[i].dot(sym_U[i]);
			sym_u[i] = Tau(q_index) - sym_S[i].dot(sym_pA[i]);
			
//			clog << "sym_D[" << i << "] = " << sym_D[i] << endl << endl;
//			clog << "sym_u[" << i << "] = " << sym_u[i] << endl << endl;

			unsigned int lambda = model.lambda[i];
			if (lambda != 0) {
				
				SymSpatialMatrix Ia =    sym_IA[i].M
				- (sym_U[i].v.transpose()
				* (sym_U[i] / sym_D[i]).v);

				SymSpatialVector pa =  sym_pA[i]
				+ Ia * sym_c[i]
				+ sym_U[i] * sym_u[i] / sym_D[i];
				
//				clog << "Ia(" << i << ").M = " << Ia.M << endl;
//				clog << "pa(" << i << ").v = " << pa.v << endl;
				
				clog << "sym_X_lambda[" << i << "].E = " << sym_X_lambda[i].E << endl;
				clog << "sym_X_lambda[" << i << "].r = " << sym_X_lambda[i].r << endl;
				
				sym_IA[lambda]
				+= sym_X_lambda[i].toMatrixTranspose()
				* (Ia * sym_X_lambda[i].toMatrix());
				
				clog << "sym_IA_2[" << lambda << "] = " << sym_IA[lambda].M << std::endl;
				
				sym_pA[lambda] += sym_X_lambda[i].applyTranspose(pa);

				clog << "sym_pA_2(" << lambda << ") = " << sym_pA[lambda] << std::endl;
				
				timer1.stop();
			}
		} else {
			std::cerr << "Error: invalid joint type " << model.mJoints[i].mJointType << " at id " << i << std::endl;
			abort();
		} /* if (model.mJoints[i].mDoFCount == 1) */
	} /* SECOND LOOP */
	
	
	sym_a[0] = spatial_gravity * Symbolic(-1.0);
	
	Symbolic QDDot("QDDot", model.dof_count);
	
	/* THIRD LOOP */
	for (i = 1; i < model.mBodies.size(); i++) {
		clog << "[LOOP 3] Iteration "<< i << "/" << model.mBodies.size() - 1 << endl << endl;
		timer1.start();
		
		unsigned int q_index = model.mJoints[i].q_index;
		unsigned int lambda = model.lambda[i];
		SymSpatialTransform X_lambda = sym_X_lambda[i];
		
//		clog << "sym_a[lambda] = " << sym_a[lambda].v << endl;

		sym_a[i] = X_lambda.apply(sym_a[lambda]) + sym_c[i];
//		clog << "a'[" << i << "] = " << sym_a[i].v << std::endl;

		if (model.mJoints[i].mDoFCount == 1
			&& model.mJoints[i].mJointType != JointTypeCustom) {
			QDDot(q_index) = (symOne/sym_D[i]) * (sym_u[i] - sym_U[i].dot(sym_a[i]));
			sym_a[i] = sym_a[i] + sym_S[i] * QDDot(q_index);
		} else {
			std::cerr << "Error: invalid joint type " << model.mJoints[i].mJointType << " at id " << i << std::endl;
			abort();
		}
		
//		clog << "QDDot(" << q_index << ") = " << QDDot(q_index) << endl;
		timer1.stop();
	} /* THIRD LOOP */
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
	} /* if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) */
	sym_X_lambda[joint_id] = sym_X_J[joint_id] * SpatialTrans2SymSpatialTrans(model.X_T[joint_id]);
}


#endif /* SymbolicForwardDynamics_hpp */
