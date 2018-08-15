//
//  HelperFunctions.hpp
//  urdf2eom
//
//  Created by Stefan on 08.08.18.
//  Copyright © 2018 Stefan Durstberger. All rights reserved.
//

#ifndef HelperFunctions_hpp
#define HelperFunctions_hpp

#include <stdio.h>
#include <chrono>  // for high_resolution_clock
#include <iostream>
#include <sstream>

class Timer {
	
private:
	int timerNumber;
	bool timerOff;
	
	std::chrono::time_point<std::__1::chrono::steady_clock, std::chrono::duration<long long, std::ratio<1LL, 1000000000LL> > > starttime;
	std::chrono::time_point<std::__1::chrono::steady_clock, std::chrono::duration<long long, std::ratio<1LL, 1000000000LL> > > endtime;
	std::chrono::time_point<std::__1::chrono::steady_clock, std::chrono::duration<long long, std::ratio<1LL, 1000000000LL> > > laptime1;
	std::chrono::time_point<std::__1::chrono::steady_clock, std::chrono::duration<long long, std::ratio<1LL, 1000000000LL> > > laptime2;
	void output (const std::ostringstream& os);
	
public:
	static int numOfTimers;
	static bool allTimersOff;
	
	Timer();
	Timer(bool timOff);
	void start (void);
	void lap (void);
	void lap (const std::string& input);
	void stop (void);
	void OFF (void);
	void ON (void);
};

#endif /* HelperFunctions_hpp */
