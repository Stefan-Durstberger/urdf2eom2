//
//  HelperFunctions.cpp
//  urdf2eom
//
//  Created by Stefan on 08.08.18.
//  Copyright © 2018 Stefan Durstberger. All rights reserved.
//

#include "HelperFunctions.hpp"

Timer::Timer() {
	timerNumber = numOfTimers;
	numOfTimers++;
	timerOff = false;
};
Timer::Timer(bool timOff) {
	timerNumber = numOfTimers;
	numOfTimers++;
	timerOff = timOff;
};
void Timer::start (void){
	starttime = std::chrono::high_resolution_clock::now();

	std::ostringstream oss;
	oss << "[TIMER " << timerNumber << "]  Started " << std::endl << std::endl;
	output (oss);
	
	laptime1 = starttime;
};
void Timer::lap (void){
	laptime2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = laptime2 - laptime1;
	
	std::ostringstream oss;
	oss << "[TIMER " << timerNumber << "]  Lap = "<< elapsed.count() << std::endl << std::endl;
	output (oss);
	
	laptime1 = laptime2;
};
void Timer::stop (void){
	endtime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = endtime - starttime;
	
	std::ostringstream oss;
	oss << "[TIMER " << timerNumber << "]  Total = "<< elapsed.count() << std::endl << std::endl;
	output (oss);
};
void Timer::output (const std::ostringstream& oss){
	if ( ! Timer::allTimersOff && ! timerOff) {
		std::cout << oss.str();
	};
};
void Timer::OFF (void) {
	timerOff = true;
};
void Timer::ON (void) {
	timerOff = false;
};
	
int Timer::numOfTimers = 0;
bool Timer::allTimersOff = false;
