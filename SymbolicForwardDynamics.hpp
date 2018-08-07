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

using namespace SymbolicRigidBodyDynamics;
using namespace Math;
using namespace std;


Symbolic symZero(0.0);
Symbolic symOne(1.0);


struct SymSpatialVector : public Symbolic {
	SymSpatialVector() :
	v{"v", 1, 6}
	{}
	SymSpatialVector (const Symbolic &vector) :
	v (vector)
	{}
	Symbolic v;
	
	SymSpatialVector operator+ (const SymSpatialVector &XT) const {
//		std::clog << "XT.v = " << XT.v.rows() << "x" << XT.v.columns() << std::endl;
//		std::clog << "v = " << v.rows() << "x" << v.columns() << std::endl;
		return SymSpatialVector (v + XT.v);
	}
}; /* struct SymSpatialVector : public Symbolic */


struct SymSpatialMatrix : public Symbolic {
	SymSpatialMatrix() :
	M{"M", 6, 6}
	{}
	SymSpatialMatrix (const Symbolic &matrix) :
	M (matrix)
	{}
	
public:
	Symbolic M;
}; /* struct SymSpatialMatrix : public Symbolic */


struct SymSpatialTransform {
	SymSpatialTransform() :
	E{"E", 3, 3},
	r{"r", 3}
	{}
	SymSpatialTransform (const Symbolic &rotation, const Symbolic &translation) :
	E (rotation),
	r (translation)
	{}
	
	
	SymSpatialTransform operator* (const SymSpatialTransform &XT) const {
		//		std::clog << "E = " << E.rows() << "x" << E.columns() << std::endl;
		//		std::clog << "XT.E = " << XT.E.rows() << "x" << XT.E.columns() << std::endl;
		//		std::clog << "r = " << r.rows() << "x" << r.columns() << std::endl;
		//		std::clog << "XT.r = " << XT.r.rows() << "x" << XT.r.columns() << std::endl;
		return SymSpatialTransform (E * XT.E, XT.r + XT.E.transpose() * r);
	}
	
	
	SymSpatialVector apply (const SymSpatialVector &v_sp) {
//		clog << "r = " << r << endl;
//		clog << "v_sp.v = " << v_sp.v << endl;
		
		Symbolic v_rxw (list<Symbolic> {
			Symbolic (v_sp.v(3) - r(1)*v_sp.v(2) + r(2)*v_sp.v(1)),
			Symbolic (v_sp.v(4) - r(2)*v_sp.v(0) + r(0)*v_sp.v(2)),
			Symbolic (v_sp.v(5) - r(0)*v_sp.v(1) + r(1)*v_sp.v(0))
		});
//		clog << "v_rxw = " << v_rxw << endl;
		
		return SymSpatialVector (Symbolic (list<Symbolic>{
			Symbolic (E(0,0) * v_sp.v(0) + E(0,1) * v_sp.v(1) + E(0,2) * v_sp.v(2)),
			Symbolic (E(1,0) * v_sp.v(0) + E(1,1) * v_sp.v(1) + E(1,2) * v_sp.v(2)),
			Symbolic (E(2,0) * v_sp.v(0) + E(2,1) * v_sp.v(1) + E(2,2) * v_sp.v(2)),
			Symbolic (E(0,0) * v_rxw(0) + E(0,1) * v_rxw(1) + E(0,2) * v_rxw(2)),
			Symbolic (E(1,0) * v_rxw(0) + E(1,1) * v_rxw(1) + E(1,2) * v_rxw(2)),
			Symbolic (E(2,0) * v_rxw(0) + E(2,1) * v_rxw(1) + E(2,2) * v_rxw(2))
		}));
	}
	

	SymSpatialMatrix toMatrixAdjoint () const {
		Symbolic _Erx =
		E * Symbolic (list<list<Symbolic>> {
			{symZero, -r(2), r(1)},
			{r(2), symZero, -r(0)},
			{-r(1), r(0), symZero}});
		
		SymSpatialMatrix result;
		for (unsigned int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j < 3; j++) {
				result.M(i,j) = E(i,j);
			}
		}
		for (unsigned int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j < 3; j++) {
				result.M(i,j+3) = - _Erx(i,j);
			}
		}
		for (unsigned int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j < 3; j++) {
				result.M(i+3,j) = symZero;
			}
		}
		for (unsigned int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j < 3; j++) {
				result.M(i+3,j+3) = E(i,j);
			}
		}
		return result;
	}
	
	Symbolic E;
	Symbolic r;
}; /* struct SymSpatialTransform */

Symbolic symVectorZero3d (list<Symbolic> {symZero, symZero, symZero});
Symbolic symVectorZero6d (list<Symbolic> {symZero, symZero, symZero, symZero, symZero, symZero});
Symbolic symMatrixZero3x3 (list<list<Symbolic>> {{symZero, symZero, symZero},
												 {symZero, symZero, symZero},
												 {symZero, symZero, symZero}});
Symbolic symMatrixZero6x6 (list<list<Symbolic>> {{symZero, symZero, symZero, symZero, symZero, symZero},
												 {symZero, symZero, symZero, symZero, symZero, symZero},
												 {symZero, symZero, symZero, symZero, symZero, symZero},
												 {symZero, symZero, symZero, symZero, symZero, symZero},
												 {symZero, symZero, symZero, symZero, symZero, symZero},
												 {symZero, symZero, symZero, symZero, symZero, symZero}});
Symbolic symMatrixIdentity6x6 (list<list<Symbolic>> {{symOne, symZero, symZero, symZero, symZero, symZero},
												 {symZero, symOne, symZero, symZero, symZero, symZero},
												 {symZero, symZero, symOne, symZero, symZero, symZero},
												 {symZero, symZero, symZero, symOne, symZero, symZero},
												 {symZero, symZero, symZero, symZero, symOne, symZero},
												 {symZero, symZero, symZero, symZero, symZero, symOne}});
SymSpatialVector symSpatialVectorZero (symVectorZero6d);
SymSpatialMatrix symSpatialMatrixZero (symMatrixZero6x6);
SymSpatialMatrix symSpatialMatrixIdentity (symMatrixIdentity6x6);


///** \brief Compact representation for Spatial Inertia. */
struct SymSpatialRigidBodyInertia {
	SymSpatialRigidBodyInertia() :
	m (symZero),
	h (symVectorZero3d),
	Ixx (symZero), Iyx(symZero), Iyy(symZero), Izx(symZero), Izy(symZero), Izz(symZero)
	{}
	SymSpatialRigidBodyInertia (
							 const Symbolic &mass, const Symbolic &com_mass,const Symbolic &inertia) :
	m (mass), h (com_mass),
	Ixx (inertia(0,0)),
	Iyx (inertia(1,0)), Iyy(inertia(1,1)),
	Izx (inertia(2,0)), Izy(inertia(2,1)), Izz(inertia(2,2))
	{ }
	SymSpatialRigidBodyInertia (const Symbolic m, const Symbolic &h,
							 const Symbolic &Ixx,
							 const Symbolic &Iyx, const Symbolic &Iyy,
							 const Symbolic &Izx, const Symbolic &Izy, const Symbolic &Izz
							 ) :
	m (m), h (h),
	Ixx (Ixx),
	Iyx (Iyx), Iyy(Iyy),
	Izx (Izx), Izy(Izy), Izz(Izz)
	{ }

	SymSpatialVector operator* (const SymSpatialVector &mv) {
		Symbolic mv_lower (list<Symbolic> {mv.v(3), mv.v(4), mv.v(5)} );
		
		clog << "mv_lower = " << mv_lower << endl;

		Symbolic res_upper (list<Symbolic> {
			Ixx * mv.v(0) + Iyx * mv.v(1) + Izx * mv.v(2),
			Iyx * mv.v(0) + Iyy * mv.v(1) + Izy * mv.v(2),
			Izx * mv.v(0) + Izy * mv.v(1) + Izz * mv.v(2)
			} + cross(h, mv_lower));
		
		clog << "res_upper = " << res_upper << endl;
		
		Symbolic res_lower = m * mv_lower - cross (h, Symbolic (list<Symbolic> {mv.v(0), mv.v(1), mv.v(2)}));
		
		clog << "res_lower = " << res_lower << endl;

		return SymSpatialVector (Symbolic (list<Symbolic> {
			res_upper(0), res_upper(1), res_upper(2),
			res_lower(0), res_lower(1), res_lower(2)
		}));
		
	};
public:
	Symbolic cross( const Symbolic &v1, const Symbolic &v2) {
		return Symbolic (list<Symbolic> {v1(1) * v2(2) - v1(2) * v2(1),
										 v1(2) * v2(0) - v1(0) * v2(2),
										 v1(0) * v2(1) - v1(1) * v2(0)});
	};
		
//	SpatialRigidBodyInertia operator+ (const SpatialRigidBodyInertia &rbi) {
//		return SpatialRigidBodyInertia (
//										m + rbi.m,
//										h + rbi.h,
//										Ixx + rbi.Ixx,
//										Iyx + rbi.Iyx, Iyy + rbi.Iyy,
//										Izx + rbi.Izx, Izy + rbi.Izy, Izz + rbi.Izz
//										);
//	}
//
//	void createFromMatrix (const SpatialMatrix &Ic) {
//		m = Ic(3,3);
//		h.set (-Ic(1,5), Ic(0,5), -Ic(0,4));
//		Ixx = Ic(0,0);
//		Iyx = Ic(1,0); Iyy = Ic(1,1);
//		Izx = Ic(2,0); Izy = Ic(2,1); Izz = Ic(2,2);
//	}
//
//	SpatialMatrix toMatrix() const {
//		SpatialMatrix result;
//		result(0,0) = Ixx; result(0,1) = Iyx; result(0,2) = Izx;
//		result(1,0) = Iyx; result(1,1) = Iyy; result(1,2) = Izy;
//		result(2,0) = Izx; result(2,1) = Izy; result(2,2) = Izz;
//
//		result.block<3,3>(0,3) = VectorCrossMatrix(h);
//		result.block<3,3>(3,0) = - VectorCrossMatrix(h);
//		result.block<3,3>(3,3) = Matrix3d::Identity(3,3) * m;
//
//		return result;
//	}
//
	void setSpatialMatrix (Symbolic &mat) const {
		mat(0,0) = Ixx; mat(0,1) = Iyx; mat(0,2) = Izx;
		mat(1,0) = Iyx; mat(1,1) = Iyy; mat(1,2) = Izy;
		mat(2,0) = Izx; mat(2,1) = Izy; mat(2,2) = Izz;

		mat(3,0) =    0.; mat(3,1) =  h(2); mat(3,2) = -h(1);
		mat(4,0) = -h(2); mat(4,1) =    0.; mat(4,2) =  h(0);
		mat(5,0) =  h(1); mat(5,1) = -h(0); mat(5,2) =    0.;

		mat(0,3) =    0.; mat(0,4) = -h(2); mat(0,5) =  h(1);
		mat(1,3) =  h(2); mat(1,4) =    0.; mat(1,5) = -h(0);
		mat(2,3) = -h(1); mat(2,4) =  h(0); mat(2,5) =    0.;

		mat(3,3) =     m; mat(3,4) =    0.; mat(3,5) =    0.;
		mat(4,3) =    0.; mat(4,4) =     m; mat(4,5) =    0.;
		mat(5,3) =    0.; mat(5,4) =    0.; mat(5,5) =     m;
	}
//
//	static SpatialRigidBodyInertia createFromMassComInertiaC (double mass, const Vector3d &com, const Matrix3d &inertia_C) {
//		SpatialRigidBodyInertia result;
//		result.m = mass;
//		result.h = com * mass;
//		Matrix3d I = inertia_C + VectorCrossMatrix (com) * VectorCrossMatrix(com).transpose() * mass;
//		result.Ixx = I(0,0);
//		result.Iyx = I(1,0);
//		result.Iyy = I(1,1);
//		result.Izx = I(2,0);
//		result.Izy = I(2,1);
//		result.Izz = I(2,2);
//		return result;
//	}

	/// Mass
	Symbolic m;
	/// Coordinates of the center of mass
	Symbolic h;
	/// Inertia expressed at the origin
	Symbolic Ixx, Iyx, Iyy, Izx, Izy, Izz;
}; /* struct SymSpatialRigidBodyInertia */

SymSpatialRigidBodyInertia symSRBIZero (0.,
							symVectorZero3d,
							symMatrixZero3x3);

static inline void SymbolicForwardDynamics ( Model &model, std::vector<SymSpatialVector> *f_ext = NULL );
static inline void jcalc ( Model &model, unsigned int joint_id );
static inline SymSpatialTransform SpatialTrans2SymSpatialTrans (SpatialTransform spat);
inline SymSpatialTransform SymXrotx (const Symbolic &xrot);
inline SymSpatialTransform SymXroty (const Symbolic &xrot);
inline SymSpatialTransform SymXrotz (const Symbolic &xrot);
inline SymSpatialVector crossm (const SymSpatialVector &v1, const SymSpatialVector &v2);
inline SymSpatialMatrix crossm (const SymSpatialVector &v);
inline SymSpatialVector crossf (const SymSpatialVector &v1, const SymSpatialVector &v2);
inline SymSpatialMatrix crossf (const SymSpatialVector &v);


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
//std::vector<Math::SpatialVector> S;

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
//std::vector<Math::SpatialVector> U;
/// \brief Temporary variable D_i (RBDA p. 130)
//Math::VectorNd d;
/// \brief Temporary variable u (RBDA p. 130)
//Math::VectorNd u;
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
std::vector<Math::SpatialMatrix> IA;

/// \brief All bodies that are attached to a body via a fixed joint.
//std::vector<FixedBody> mFixedBodies;


static inline void SymbolicForwardDynamics ( Model &model, std::vector<SymSpatialVector> *f_ext ) {
	
	for (unsigned int i = 0; i < model.dof_count + 1; i++) {
		sym_X_J.push_back 		(SymSpatialTransform());
		sym_X_lambda.push_back	(SymSpatialTransform());
		sym_X_base.push_back 	(SymSpatialTransform());
		sym_v_J.push_back 		(symSpatialVectorZero);
		sym_v.push_back 		(symSpatialVectorZero);
		sym_c_J.push_back 		(symSpatialVectorZero);
		sym_c.push_back 		(symSpatialVectorZero);
		sym_IA.push_back		(symSpatialMatrixIdentity);
		//sym_Ic.push_back		(symSRBIZero);
		sym_I.push_back			(symSRBIZero);
		sym_pA.push_back		(symSpatialVectorZero);
	}
	
	//clog << symSpatialVectorZero.v << endl;
	
	SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);
	
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
		unsigned int lambda = model.lambda[i];
		
		jcalc (model, i);
		
		if (lambda != 0) {
			sym_X_base[i] = sym_X_lambda[i] * SpatialTrans2SymSpatialTrans (model.X_base[lambda]);
		}
		else
			sym_X_base[i] = sym_X_lambda[i];
		
//		clog << "sym_X_lambda[" << i << "].E = " << sym_X_lambda[i].E << endl;
//		clog << "sym_X_lambda[" << i << "].r = " << sym_X_lambda[i].r << endl;
//		clog << "sym_v[" << lambda << "].v = " << sym_v[lambda].v << endl;
//		clog << "sym_v_J[" << i << "].v = " << sym_v_J[i].v << endl;
		
		sym_v[i] = sym_X_lambda[i].apply( sym_v[lambda].v) + sym_v_J[i];
		
//		clog << "sym_v[" << i << "].v = " << sym_v[i].v << endl;
		
		sym_c[i] = sym_c_J[i] + crossm(sym_v[i],sym_v_J[i]);
		
		sym_I[i].setSpatialMatrix (sym_IA[i].M);

		sym_pA[i] = crossf(sym_v[i].v,sym_I[i] * sym_v[i]);

		if (f_ext != NULL && (*f_ext)[i] != symSpatialVectorZero) {
			clog << "External force (" << i << ") = " << sym_X_base[i].toMatrixAdjoint() * (*f_ext)[i] << std::endl;
			sym_pA[i] -= sym_X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
	}
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
		sym_v_J[joint_id][0] = qdot(model.mJoints[joint_id].q_index);
//		std::clog << "[jcalc]: sym_X_J = " << std::endl;
//		std::clog << sym_X_J[joint_id].E << std::endl;
//		std::clog << sym_X_J[joint_id].r << std::endl;
	} else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
		sym_X_J[joint_id] = SymXroty (q(model.mJoints[joint_id].q_index));
		sym_v_J[joint_id][1] = qdot(model.mJoints[joint_id].q_index);
	} else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
		sym_X_J[joint_id] = SymXrotz (q(model.mJoints[joint_id].q_index));
		sym_v_J[joint_id][3] = qdot(model.mJoints[joint_id].q_index);
	/* Joint Helical, Sherical, EulerZYX, EulerXYZ, EulerYXZ, TranlationalXYZ, Custom not available */
	} else {
		std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
		abort();
	}
	sym_X_lambda[joint_id] = sym_X_J[joint_id] * SpatialTrans2SymSpatialTrans(model.X_T[joint_id]);
}


static inline SymSpatialTransform SpatialTrans2SymSpatialTrans (SpatialTransform spat) {
	Symbolic E("E", static_cast<unsigned int>(spat.E.rows()), static_cast<unsigned int>(spat.E.cols()));
	Symbolic r("r", static_cast<unsigned int>(spat.r.size()));
	//	clog << "E = " << E << endl;
	//	clog << "r = " << r << endl;
	//	clog << "spat.E = " << spat.E << endl;
	//	clog << "spat.r = " << spat.r << endl;
	for (unsigned int i = 0; i < spat.E.rows(); i++) {
		for (unsigned int j = 0; j < spat.E.cols(); j++) {
			E(i,j) = spat.E(i, j);
		}
	}
	for (unsigned int i = 0; i < spat.r.size(); i++) {
		r(i) = spat.r[i];
	}
	//	clog << "E = " << E.rows() << E.columns() << endl;
	//	clog << "r = " << r.rows() << r.columns() << endl;
	return SymSpatialTransform ( E, r );
}



inline SymSpatialTransform SymXrotx (const Symbolic &xrot) {
	Symbolic s, c;
	Symbolic one(1.0);
	Symbolic zero(0.0);
	s = sin (xrot);
	c = cos (xrot);
	
	Symbolic M = ( (one,   zero,  zero),
				  (zero,  c,     s ),
				  (zero,  -s,    c ));
	Symbolic v = (  zero,  zero,  zero);
	//	SymSpatialTransform tmp = SymSpatialTransform(M, v);
	//	std::cout << tmp.E << std::endl;
	//	std::cout << tmp.r << std::endl;
	
	return SymSpatialTransform ( M, v.transpose() );
}



inline SymSpatialTransform SymXroty (const Symbolic &xrot) {
	Symbolic s, c;
	Symbolic symOne(1.0);
	Symbolic symZero(0.0);
	s = sin (xrot);
	c = cos (xrot);
	
	Symbolic M = ( (c,         symZero,  -s),
				  (symZero,  symOne,   symZero ),
				  (s,         symZero,  c ));
	Symbolic v = (  symZero,  symZero,  symZero);
	//	SymSpatialTransform tmp = SymSpatialTransform(M, v);
	//	std::cout << tmp.E << std::endl;
	//	std::cout << tmp.r << std::endl;
	
	return SymSpatialTransform ( M, v.transpose() );
}



inline SymSpatialTransform SymXrotz (const Symbolic &xrot) {
	Symbolic s, c;
	Symbolic one(1.0);
	Symbolic zero(0.0);
	s = sin (xrot);
	c = cos (xrot);
	
	Symbolic M = ( (c,     s,     zero),
				  (-s,    c,     zero ),
				  (zero,  zero,  one ));
	Symbolic v = (  zero,  zero,  zero);
	//	SymSpatialTransform tmp = SymSpatialTransform(M, v);
	//	std::cout << tmp.E << std::endl;
	//	std::cout << tmp.r << std::endl;
	
	return SymSpatialTransform ( M, v.transpose() );
}


inline SymSpatialMatrix crossm (const SymSpatialVector &v) {
	return SymSpatialMatrix (Symbolic (list<list<Symbolic>>{
		{symZero,  -v.v(2),  v.v(1),         symZero,          symZero,         symZero},
		{v.v(2),          symZero, -v.v(0),         symZero,          symZero,         symZero},
		{-v.v(1),   v.v(0),         symZero,         symZero,          symZero,         symZero},
		{symZero,  -v.v(5),  v.v(4),         symZero,  -v.v(2),  v.v(1)},
		{v.v(5),          symZero, -v.v(3),  v.v(2),          symZero, -v.v(0)},
		{-v.v(4),   v.v(3),         symZero, -v.v(1),   v.v(0),         symZero}}));
};


inline SymSpatialVector crossm (const SymSpatialVector &v1, const SymSpatialVector &v2) {
	return SymSpatialVector (Symbolic (list<Symbolic>{
		Symbolic (-v1.v(2) * v2.v(1) + v1.v(1) * v2.v(2)),
		Symbolic (v1.v(2) * v2.v(0) - v1.v(0) * v2.v(2)),
		Symbolic (-v1.v(1) * v2.v(0) + v1.v(0) * v2.v(1)),
		Symbolic (-v1.v(5) * v2.v(1) + v1.v(4) * v2.v(2) - v1.v(2) * v2.v(4) + v1.v(1) * v2.v(5)),
		Symbolic (v1.v(5) * v2.v(0) - v1.v(3) * v2.v(2) + v1.v(2) * v2.v(3) - v1.v(0) * v2.v(5)),
		Symbolic (-v1.v(4) * v2.v(0) + v1.v(3) * v2.v(1) - v1.v(1) * v2.v(3) + v1.v(0) * v2.v(4))
	}));
}


inline SymSpatialMatrix crossf (const SymSpatialVector &v) {
	return SymSpatialMatrix (Symbolic (list<list<Symbolic>>{
		{symZero,  -v.v(2),  v.v(1),         symZero,  -v.v(5),  v.v(4)},
		{v.v(2),          symZero, -v.v(0),  v.v(5),          symZero, -v.v(3)},
		{-v.v(1),   v.v(0),         symZero, -v.v(4),   v.v(3),         symZero},
		{symZero,          symZero,         symZero,         symZero,  -v.v(2),  v.v(1)},
		{symZero,          symZero,         symZero,  v.v(2),          symZero, -v.v(0)},
		{symZero,          symZero,         symZero, -v.v(1),   v.v(0),         symZero}}));
}


inline SymSpatialVector crossf (const SymSpatialVector &v1, const SymSpatialVector &v2) {
	return SymSpatialVector (Symbolic (list<Symbolic>{
		Symbolic (-v1.v(2) * v2.v(1) + v1.v(1) * v2.v(2) - v1.v(5) * v2.v(4) + v1.v(4) * v2.v(5)),
		Symbolic (v1.v(2) * v2.v(0) - v1.v(0) * v2.v(2) + v1.v(5) * v2.v(3) - v1.v(3) * v2.v(5)),
		Symbolic (-v1.v(1) * v2.v(0) + v1.v(0) * v2.v(1) - v1.v(4) * v2.v(3) + v1.v(3) * v2.v(4)),
		Symbolic (- v1.v(2) * v2.v(4) + v1.v(1) * v2.v(5)),
		Symbolic (+ v1.v(2) * v2.v(3) - v1.v(0) * v2.v(5)),
		Symbolic ( - v1.v(1) * v2.v(3) + v1.v(0) * v2.v(4))}));
}

#endif /* SymbolicForwardDynamics_hpp */
