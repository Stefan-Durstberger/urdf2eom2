/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <stdio.h>
#include <iostream>
#include <chrono>  // for high_resolution_clock

#include <srbdl.h>
#include <srbdl/addons/urdfreader/urdfreader.h>

#include "symbolicc++/symbolicc++.h"


#ifndef SRBDL_BUILD_ADDON_URDFREADER
	#error "Error: SRBDL addon URDFReader not enabled."
#endif

using namespace SymbolicRigidBodyDynamics;
using namespace SymbolicRigidBodyDynamics::Math;

int main (int argc, char* argv[]) {
	
//	// Record start time
//	auto start = std::chrono::high_resolution_clock::now();
	
	// Symbolic Forward Dynamics
	Model* model = new Model();

	if ( argc <= 1 ) {
		std::cerr << "No input argument" << std::endl;
		abort();
	}
	
	if (!Addons::URDFReadFromFile (argv[1], model, false)) {
		std::cerr << "Error loading model ./samplemodel.urdf" << std::endl;
		abort();
	}

	std::cout << "Degree of freedom: " << model->dof_count << std::endl;
	
	
	
	delete model;
	
	
//	// Record end time
//	auto finish = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> elapsed = finish - start;
//	std::cout << "Time elapsed = "<< elapsed.count() << std::endl;
	
	return 0;
}

