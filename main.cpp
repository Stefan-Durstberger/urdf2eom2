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

#include "SymbolicForwardDynamics.hpp"


#ifndef SRBDL_BUILD_ADDON_URDFREADER
	#error "Error: SRBDL addon URDFReader not enabled."
#endif

#define SRBDL_MEASURE_TIME


using namespace SymbolicRigidBodyDynamics;
using namespace SymbolicRigidBodyDynamics::Math;

int main (int argc, char* argv[]) {
	
	#ifdef SRBDL_MEASURE_TIME
	// Record start time
	auto start = std::chrono::high_resolution_clock::now();
	#endif
	
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
	
	
	
	
	
	
	// Symbolic Forward Dynamics
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
	
	delete model1;
	delete model2;
	
	#ifdef SRBDL_MEASURE_TIME
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Time elapsed = "<< elapsed.count() << std::endl;
	#endif
	
	return 0;
}

