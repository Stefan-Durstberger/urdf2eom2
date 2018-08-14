/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <stdio.h>
#include <iostream>

#include <srbdl.h>
#include <srbdl/addons/urdfreader/urdfreader.h>
#include "symbolicc++/symbolicc++.h"

#include "SymbolicForwardDynamics.hpp"

#include "HelperFunctions.hpp"


#ifndef SRBDL_BUILD_ADDON_URDFREADER
	#error "Error: SRBDL addon URDFReader not enabled."
#endif

using namespace SymbolicRigidBodyDynamics;
using namespace SymbolicRigidBodyDynamics::Math;

int main (int argc, char* argv[]) {
	
	// Record start time
	Timer timer; timer.start();
	
	// Symbolic Forward Dynamics
	Model* model1 = new Model();

	if ( argc <= 1 ) {
		std::cerr << "No input argument" << std::endl;
		abort();
	}
	
	if (!Addons::URDFReadFromFile (argv[1], model1, false)) {
		std::cerr << "Error loading model" << std::endl;
		abort();
	}

	//std::cout << "Degree of freedom: " << model->dof_count << std::endl;
	
	SymbolicForwardDynamics(*model1);
	
	timer.lap();
	
	// Symbolic Forward Dynamics
	cout << "[EVALUATED MODEL]" << endl;
	Model* model2 = new Model();
	
	if ( argc <= 1 ) {
		std::cerr << "No input argument" << std::endl;
		abort();
	}
	
	if (!Addons::URDFReadFromFile (argv[1], model2, false)) {
		std::cerr << "Error loading model" << std::endl;
		abort();
	}
	
	VectorNd Q = VectorNd::Ones (model2->dof_count);
	VectorNd QDot = VectorNd::Ones (model2->dof_count);
	VectorNd Tau = VectorNd::Ones (model2->dof_count);
	VectorNd QDDot = VectorNd::Ones (model2->dof_count);
	
	ForwardDynamics (*model2, Q, QDot, Tau, QDDot);
	
	cout << "QDDot = " << endl << QDDot << endl << endl;
	
	delete model1;
	delete model2;
	
	// Record end time
	timer.stop();
	
	return 0;
}

