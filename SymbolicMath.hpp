//
//  SymbolicMathUtil.hpp
//  urdf2eom
//
//  Created by Stefan on 08.08.18.
//  Copyright © 2018 Stefan Durstberger. All rights reserved.
//

#ifndef SymbolicMathUtil_hpp
#define SymbolicMathUtil_hpp

#include <stdio.h>
#include "SymbolicDataTypes.hpp"
#include "srbdl/Model.h"

namespace SymbolicMath {

using namespace SymbolicDataTypes;

static inline SymSpatialTransform SpatialTrans2SymSpatialTrans (SymbolicRigidBodyDynamics::Math::SpatialTransform spat);
inline SymSpatialTransform SymXrotx (const Symbolic &xrot);
inline SymSpatialTransform SymXroty (const Symbolic &xrot);
inline SymSpatialTransform SymXrotz (const Symbolic &xrot);
inline SymSpatialVector crossm (const SymSpatialVector &v1, const SymSpatialVector &v2);
inline SymSpatialMatrix crossm (const SymSpatialVector &v);
inline SymSpatialVector crossf (const SymSpatialVector &v1, const SymSpatialVector &v2);
inline SymSpatialMatrix crossf (const SymSpatialVector &v);
	
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
	return SymSpatialVector ( list<Symbolic>{
		(-v1.v(2) * v2.v(1) + v1.v(1) * v2.v(2)),
		(v1.v(2) * v2.v(0) - v1.v(0) * v2.v(2)),
		(-v1.v(1) * v2.v(0) + v1.v(0) * v2.v(1)),
		(-v1.v(5) * v2.v(1) + v1.v(4) * v2.v(2) - v1.v(2) * v2.v(4) + v1.v(1) * v2.v(5)),
		(v1.v(5) * v2.v(0) - v1.v(3) * v2.v(2) + v1.v(2) * v2.v(3) - v1.v(0) * v2.v(5)),
		(-v1.v(4) * v2.v(0) + v1.v(3) * v2.v(1) - v1.v(1) * v2.v(3) + v1.v(0) * v2.v(4))
	});
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
	
}

#endif /* SymbolicMathUtil_hpp */
