//
//  SymbolicMath.hpp
//  urdf2eom
//
//  Created by Stefan on 08.08.18.
//  Copyright © 2018 Stefan Durstberger. All rights reserved.
//

#ifndef SymbolicMath_hpp
#define SymbolicMath_hpp

#include <stdio.h>
#include "symbolicc++/symbolicc++.h"

namespace SymbolicDataTypes {
	
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
//	SymSpatialVector operator- (const SymSpatialVector &XT) const {
//		//		std::clog << "XT.v = " << XT.v.rows() << "x" << XT.v.columns() << std::endl;
//		//		std::clog << "v = " << v.rows() << "x" << v.columns() << std::endl;
//		return SymSpatialVector (v - XT.v);
//	}
	SymSpatialVector operator* (const Symbolic &XT) const {
		//		std::clog << "XT.v = " << XT.v.rows() << "x" << XT.v.columns() << std::endl;
		//		std::clog << "v = " << v.rows() << "x" << v.columns() << std::endl;
		return SymSpatialVector (v * XT);
	}
	SymSpatialVector operator* (const SymSpatialVector &XT) const {
		//		std::clog << "XT.v = " << XT.v.rows() << "x" << XT.v.columns() << std::endl;
		//		std::clog << "v = " << v.rows() << "x" << v.columns() << std::endl;
		return Symbolic (v * XT.v);
	}
	SymSpatialVector operator/ (const Symbolic &XT) const {
		//		std::clog << "XT.v = " << XT.v.rows() << "x" << XT.v.columns() << std::endl;
		//		std::clog << "v = " << v.rows() << "x" << v.columns() << std::endl;
		return SymSpatialVector (v / XT);
	}
	
	Symbolic dot (const SymSpatialVector &XT) const {
		Symbolic tmp = (v(0) * XT.v(0) + v(1) * XT.v(1) + v(2) * XT.v(2) + v(3) * XT.v(3) + v(4) * XT.v(4) + v(5) * XT.v(5));
		return tmp;
	}
	
	SymSpatialVector transpose () const {
		return SymSpatialVector (v.transpose());
	}
}; /* struct SymSpatialVector : public Symbolic */


struct SymSpatialMatrix : Symbolic {
	SymSpatialMatrix() :
	M{"M", 6, 6}
	{}
	SymSpatialMatrix (const Symbolic &matrix) :
	M (matrix)
	{}
	Symbolic M;
	
	SymSpatialVector operator* (const SymSpatialVector &XT) const {
		//		std::clog << "E = " << E.rows() << "x" << E.columns() << std::endl;
		//		std::clog << "XT.E = " << XT.E.rows() << "x" << XT.E.columns() << std::endl;
		//		std::clog << "r = " << r.rows() << "x" << r.columns() << std::endl;
		//		std::clog << "XT.r = " << XT.r.rows() << "x" << XT.r.columns() << std::endl;
		return SymSpatialVector (M * XT.v.transpose()).v.transpose();
	};
	SymSpatialMatrix operator* (const SymSpatialMatrix &XT) const {
		//		std::clog << "XT.v = " << XT.v.rows() << "x" << XT.v.columns() << std::endl;
		//		std::clog << "v = " << v.rows() << "x" << v.columns() << std::endl;
		return SymSpatialMatrix (M * XT.M);
	}
}; /* struct SymSpatialMatrix : public Symbolic */


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
	
	/// Mass
	Symbolic m;
	/// Coordinates of the center of mass
	Symbolic h;
	/// Inertia expressed at the origin
	Symbolic Ixx, Iyx, Iyy, Izx, Izy, Izz;
	
	SymSpatialVector operator* (const SymSpatialVector &mv) {
		Symbolic mv_lower (list<Symbolic> {mv.v(3), mv.v(4), mv.v(5)} );
		//clog << "mv_lower = " << mv_lower << endl;
		Symbolic res_upper (list<Symbolic> {
			Ixx * mv.v(0) + Iyx * mv.v(1) + Izx * mv.v(2),
			Iyx * mv.v(0) + Iyy * mv.v(1) + Izy * mv.v(2),
			Izx * mv.v(0) + Izy * mv.v(1) + Izz * mv.v(2)
		} + cross(h, mv_lower));
		//clog << "res_upper = " << res_upper << endl;
		Symbolic res_lower = m * mv_lower - cross (h, Symbolic (list<Symbolic> {mv.v(0), mv.v(1), mv.v(2)}));
		//clog << "res_lower = " << res_lower << endl;
		return SymSpatialVector (Symbolic (list<Symbolic> {
			res_upper(0), res_upper(1), res_upper(2),
			res_lower(0), res_lower(1), res_lower(2)
		}));
		
	};
	
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
	void setSpatialMatrix (SymSpatialMatrix &mat) const {
		mat.M(0,0) = Ixx; mat.M(0,1) = Iyx; mat.M(0,2) = Izx;
		mat.M(1,0) = Iyx; mat.M(1,1) = Iyy; mat.M(1,2) = Izy;
		mat.M(2,0) = Izx; mat.M(2,1) = Izy; mat.M(2,2) = Izz;
		
		mat.M(3,0) =    0.; mat.M(3,1) =  h(2); mat.M(3,2) = -h(1);
		mat.M(4,0) = -h(2); mat.M(4,1) =    0.; mat.M(4,2) =  h(0);
		mat.M(5,0) =  h(1); mat.M(5,1) = -h(0); mat.M(5,2) =    0.;
		
		mat.M(0,3) =    0.; mat.M(0,4) = -h(2); mat.M(0,5) =  h(1);
		mat.M(1,3) =  h(2); mat.M(1,4) =    0.; mat.M(1,5) = -h(0);
		mat.M(2,3) = -h(1); mat.M(2,4) =  h(0); mat.M(2,5) =    0.;
		
		mat.M(3,3) =     m; mat.M(3,4) =    0.; mat.M(3,5) =    0.;
		mat.M(4,3) =    0.; mat.M(4,4) =     m; mat.M(4,5) =    0.;
		mat.M(5,3) =    0.; mat.M(5,4) =    0.; mat.M(5,5) =     m;
	}
	
	static SymSpatialRigidBodyInertia createFromMassComInertiaC (const Symbolic &mass, const Symbolic &com,const Symbolic &inertia_C) {
		SymSpatialRigidBodyInertia result;
		result.m (mass);
		result.h (com * mass);
		Symbolic I = inertia_C + SymVectorCrossMatrix (com) * SymVectorCrossMatrix(com).transpose() * mass;
		result.Ixx = I(0,0);
		result.Iyx = I(1,0);
		result.Iyy = I(1,1);
		result.Izx = I(2,0);
		result.Izy = I(2,1);
		result.Izz = I(2,2);
		return result;
	}
	
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
	
	static Symbolic SymVectorCrossMatrix (const Symbolic &vector) {
		return Symbolic (list<list<Symbolic>> {
			{symZero, -vector(2), vector(1)},
			{vector(2), symZero, -vector(0)},
			{-vector(1), vector(0), symZero}});
	}
}; /* struct SymSpatialRigidBodyInertia */

SymSpatialRigidBodyInertia symSRBIZero (0.,
										symVectorZero3d,
										symMatrixZero3x3);

struct SymSpatialTransform {
	SymSpatialTransform() :
	E{"E", 3, 3},
	r{"r", 3}
	{}
	SymSpatialTransform (const Symbolic &rotation, const Symbolic &translation) :
	E (rotation),
	r (translation)
	{}
	
	Symbolic E;
	Symbolic r;
	
	/** Same as X * v.
	 *
	 * \returns (E * w, - E * rxw + E * v)
	 */
	SymSpatialVector apply (const SymSpatialVector &v_sp) {
//				clog << "r = " << r << endl;
//				clog << "v_sp.v = " << v_sp.v << endl;
		Symbolic v_rxw (list<Symbolic> {
			(v_sp.v(3) - r(1)*v_sp.v(2) + r(2)*v_sp.v(1)),
			(v_sp.v(4) - r(2)*v_sp.v(0) + r(0)*v_sp.v(2)),
			(v_sp.v(5) - r(0)*v_sp.v(1) + r(1)*v_sp.v(0))
		});
		//		clog << "v_rxw = " << v_rxw << endl;
		return SymSpatialVector (Symbolic (list<Symbolic>{
			(E(0,0) * v_sp.v(0) + E(0,1) * v_sp.v(1) + E(0,2) * v_sp.v(2)),
			(E(1,0) * v_sp.v(0) + E(1,1) * v_sp.v(1) + E(1,2) * v_sp.v(2)),
			(E(2,0) * v_sp.v(0) + E(2,1) * v_sp.v(1) + E(2,2) * v_sp.v(2)),
			(E(0,0) * v_rxw(0) + E(0,1) * v_rxw(1) + E(0,2) * v_rxw(2)),
			(E(1,0) * v_rxw(0) + E(1,1) * v_rxw(1) + E(1,2) * v_rxw(2)),
			(E(2,0) * v_rxw(0) + E(2,1) * v_rxw(1) + E(2,2) * v_rxw(2))
		}));
	};
	
	/** Same as X^T * f.
	 *
	 * \returns (E^T * n + rx * E^T * f, E^T * f)
	 */
	SymSpatialVector applyTranspose (const SymSpatialVector &f_sp) {
		Symbolic E_T_f (list<Symbolic> {
						E(0,0) * f_sp.v(3) + E(1,0) * f_sp.v(4) + E(2,0) * f_sp.v(5),
						E(0,1) * f_sp.v(3) + E(1,1) * f_sp.v(4) + E(2,1) * f_sp.v(5),
						E(0,2) * f_sp.v(3) + E(1,2) * f_sp.v(4) + E(2,2) * f_sp.v(5)});
		//		clog << "E_T_f = " << E_T_f << endl;
		return SymSpatialVector (Symbolic( list<Symbolic> {
							E(0,0) * f_sp.v(0) + E(1,0) * f_sp.v(1) + E(2,0) * f_sp.v(2) - r(2) * E_T_f(1) + r(1) * E_T_f(2),
							E(0,1) * f_sp.v(0) + E(1,1) * f_sp.v(1) + E(2,1) * f_sp.v(2) + r(2) * E_T_f(0) - r(0) * E_T_f(2),
							E(0,2) * f_sp.v(0) + E(1,2) * f_sp.v(1) + E(2,2) * f_sp.v(2) - r(1) * E_T_f(0) + r(0) * E_T_f(1),
							E_T_f (0),
							E_T_f (1),
							E_T_f (2)}));
	};
	
//	/** Same as X^* I X^{-1}
//	 */
//	SymSpatialRigidBodyInertia apply (const SymSpatialRigidBodyInertia &rbi) {
//		return SymSpatialRigidBodyInertia (
//										rbi.m,
//										E * (rbi.h - rbi.m * r),
//										E * (Symbolic (list<list<Symbolic>> {
//											{rbi.Ixx, rbi.Iyx, rbi.Izx},
//											{rbi.Iyx, rbi.Iyy, rbi.Izy},
//											{rbi.Izx, rbi.Izy, rbi.Izz}})
//										+ SymVectorCrossMatrix (r) * SymVectorCrossMatrix (rbi.h)
//										+ (SymVectorCrossMatrix(rbi.h - rbi.m * r) * SymVectorCrossMatrix (r))
//										 )
//										* E.transpose()
//										);
//	};
	
	/** Same as X^T I X
	 */
	//	SymSpatialRigidBodyInertia applyTranspose (const SymSpatialRigidBodyInertia &rbi) {
	//		Symbolic E_T_mr = E.transpose() * rbi.h + rbi.m * r;
	//		return SymSpatialRigidBodyInertia (
	//										   rbi.m,
	//										   E_T_mr,
	//										   E.transpose() *
	//										   Symbolic (list<list<Symbolic>> {
	//											{rbi.Ixx, rbi.Iyx, rbi.Izx},
	//											{rbi.Iyx, rbi.Iyy, rbi.Izy},
	//											{rbi.Izx, rbi.Izy, rbi.Izz}})
	//										   ) * E
	//										   - SymVectorCrossMatrix(r) * SymVectorCrossMatrix (E.transpose() * rbi.h)
	//										   - SymVectorCrossMatrix (E_T_mr) * SymVectorCrossMatrix (r);
	//	};
	
	//	SpatialVector applyAdjoint (const SpatialVector &f_sp) {
	//		Vector3d En_rxf = E * (Vector3d (f_sp[0], f_sp[1], f_sp[2]) - r.cross(Vector3d (f_sp[3], f_sp[4], f_sp[5])));
	//		//		Vector3d En_rxf = E * (Vector3d (f_sp[0], f_sp[1], f_sp[2]) - r.cross(Eigen::Map<Vector3d> (&(f_sp[3]))));
	//
	//		return SpatialVector (
	//							  En_rxf[0],
	//							  En_rxf[1],
	//							  En_rxf[2],
	//							  E(0,0) * f_sp[3] + E(0,1) * f_sp[4] + E(0,2) * f_sp[5],
	//							  E(1,0) * f_sp[3] + E(1,1) * f_sp[4] + E(1,2) * f_sp[5],
	//							  E(2,0) * f_sp[3] + E(2,1) * f_sp[4] + E(2,2) * f_sp[5]
	//							  );
	//	};
	
	SymSpatialMatrix toMatrix () const {
		Symbolic _Erx =
		E * Symbolic (list<list<Symbolic>> {
			{symZero, -r(2), r(1)},
			{r(2), symZero, -r(0)},
			{-r(1), r(0), symZero}});
		SymSpatialMatrix result;
		// left top = E
		for (unsigned int i = 0; i < 3; i++) //Column
			for (unsigned int j = 0; j < 3; j++) //Row
				result.M(i,j) = E(i,j);
		// left bottom = zero
		for (unsigned int i = 0; i < 3; i++) //Column
			for (unsigned int j = 0; j < 3; j++) //Row
				result.M(i,j+3) = symZero;
		// right top = -_Erx
		for (unsigned int i = 0; i < 3; i++) //Column
			for (unsigned int j = 0; j < 3; j++) //Row
				result.M(i+3,j) = -_Erx(i,j);
		// right bottom = E
		for (unsigned int i = 0; i < 3; i++) //Column
			for (unsigned int j = 0; j < 3; j++) //Row
				result.M(i+3,j+3) = E(i,j);
		return result;
	}
	
	SymSpatialMatrix toMatrixAdjoint () const {
		Symbolic _Erx =
		E * Symbolic (list<list<Symbolic>> {
			{symZero, -r(2), r(1)},
			{r(2), symZero, -r(0)},
			{-r(1), r(0), symZero}});
		SymSpatialMatrix result;
		// left top = E
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i,j) = E(i,j);
		// left bottom = -_Erx
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i,j+3) = - _Erx(i,j);
		// right top = zero
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i+3,j) = symZero;
		// right bottom = E
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i+3,j+3) = E(i,j);
		return result;
	};
	
	SymSpatialMatrix toMatrixTranspose () const {
		Symbolic _Erx =
		E * Symbolic (list<list<Symbolic>> {
			{symZero, -r(2), r(1)},
			{r(2), symZero, -r(0)},
			{-r(1), r(0), symZero}});
//		cout << "_Erx = " << _Erx << endl;
		SymSpatialMatrix result;
		// left top = E.transpose()
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i,j) = E(j,i);
		// left bottom = -_Erx.transpose
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i,j+3) = -_Erx(j,i);
		// right top = zero
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i+3,j) = symZero;
		// right bottom = E.transpose
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				result.M(i+3,j+3) = E(j,i);
		return result;
	}
	
	//	SpatialTransform inverse() const {
	//		return SpatialTransform (
	//								 E.transpose(),
	//								 - E * r
	//								 );
	//	}
	
	SymSpatialTransform operator* (const SymSpatialTransform &XT) const {
		//		std::clog << "E = " << E.rows() << "x" << E.columns() << std::endl;
		//		std::clog << "XT.E = " << XT.E.rows() << "x" << XT.E.columns() << std::endl;
		//		std::clog << "r = " << r.rows() << "x" << r.columns() << std::endl;
		//		std::clog << "XT.r = " << XT.r.rows() << "x" << XT.r.columns() << std::endl;
		return SymSpatialTransform (E * XT.E, XT.r + XT.E.transpose() * r);
	};
	
	//	void operator*= (const SpatialTransform &XT) {
	//		r = XT.r + XT.E.transpose() * r;
	//		E *= XT.E;
	//	}

	static Symbolic SymVectorCrossMatrix (const Symbolic &vector) {
		return Symbolic (list<list<Symbolic>> {
			{symZero, -vector(2), vector(1)},
			{vector(2), symZero, -vector(0)},
			{-vector(1), vector(0), symZero}});
	}
}; /* struct SymSpatialTransform */
	
} /* namespace SymbolicDataTypes */
#endif /* SymbolicMath_hpp */
