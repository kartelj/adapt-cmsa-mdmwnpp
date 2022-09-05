/***************************************************************************
                          Timer.cc  -  description
                             -------------------
    begin                : Fri Nov 10 2000
    copyright            : (C) 2000 by Christian Blum
    email                : cblum@ulb.ac.be
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "Timer.h"
#ifdef _WIN32
#include <chrono>
#endif
/*
 *  The virtual time of day and the real time of day are calculated and
 *  stored for future use.  The future use consists of subtracting these
 *  values from similar values obtained at a later time to allow the user
 *  to get the amount of time used by the algorithm.
 */
Timer::Timer(void) {
#ifdef _WIN32
    real_time = chrono::steady_clock::now();
#else
  getrusage( RUSAGE_SELF, &res );
  virtual_time = (double) res.ru_utime.tv_sec +
    (double) res.ru_stime.tv_sec +
    (double) res.ru_utime.tv_usec * 1.0E-6 +
    (double) res.ru_stime.tv_usec * 1.0E-6;
  
  gettimeofday( &tp, NULL );
  real_time =    (double) tp.tv_sec + (double) tp.tv_usec * 1.0E-6;
#endif
}

/*
 *  Stop the stopwatch and return the time used in seconds (either
 *  REAL or VIRTUAL time, depending on ``type'').
 */
double Timer::elapsed_time(const TYPE& type) {
#ifdef _WIN32
    return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - real_time).count()/1000.0;
#else
  if (type == REAL) {
    gettimeofday( &tp, NULL );
    return( (double) tp.tv_sec + (double) tp.tv_usec * 1.0E-6 - real_time );
  }
  else {
    getrusage( RUSAGE_SELF, &res );
    return( (double) res.ru_utime.tv_sec +
	    (double) res.ru_stime.tv_sec +
	    (double) res.ru_utime.tv_usec * 1.0E-6 +
	    (double) res.ru_stime.tv_usec * 1.0E-6
	    - virtual_time );
  }
#endif
}
