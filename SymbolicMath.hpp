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
#include "srbdl/Model.h"

#include "symbolicc++/symbolicc++.h"

//using namespace SymbolicRigidBodyDynamics;
//using namespace Math;
using namespace std;

namespace SymbolicMath {

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


static inline SymSpatialTransform SpatialTrans2SymSpatialTrans (SymbolicRigidBodyDynamics::Math::SpatialTransform spat);
inline SymSpatialTransform SymXrotx (const Symbolic &xrot);
inline SymSpatialTransform SymXroty (const Symbolic &xrot);
inline SymSpatialTransform SymXrotz (const Symbolic &xrot);
inline SymSpatialVector crossm (const SymSpatialVector &v1, const SymSpatialVector &v2);
inline SymSpatialMatrix crossm (const SymSpatialVector &v);
inline SymSpatialVector crossf (const SymSpatialVector &v1, const SymSpatialVector &v2);
inline SymSpatialMatrix crossf (const SymSpatialVector &v);
inline Symbolic SymVectorCrossMatrix (const Symbolic &vector);


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





static inline SymSpatialTransform SpatialTrans2SymSpatialTrans (SymbolicRigidBodyDynamics::Math::SpatialTransform spat) {
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
	s = sin (xrot);
	c = cos (xrot);
	
	Symbolic M = ( (symOne,   symZero,  symZero),
				  (symZero,  c,     s ),
				  (symZero,  -s,    c ));
	Symbolic v = (  symZero,  symZero,  symZero);
	//	SymSpatialTransform tmp = SymSpatialTransform(M, v);
	//	std::cout << tmp.E << std::endl;
	//	std::cout << tmp.r << std::endl;
	
	return SymSpatialTransform ( M, v.transpose() );
}



inline SymSpatialTransform SymXroty (const Symbolic &xrot) {
	Symbolic s, c;
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
	s = sin (xrot);
	c = cos (xrot);
	
	Symbolic M = ( (c,     s,     symZero),
				  (-s,    c,     symZero ),
				  (symZero,  symZero,  symOne ));
	Symbolic v = (  symZero,  symZero,  symZero);
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


inline Symbolic SymVectorCrossMatrix (const Symbolic &vector) {
	return Symbolic (list<list<Symbolic>> {
		{symZero, -vector(2), vector(1)},
		{vector(2), symZero, -vector(0)},
		{-vector(1), vector(0), symZero}});
}

}
#endif /* SymbolicMath_hpp */
