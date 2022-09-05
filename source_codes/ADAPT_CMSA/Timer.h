/***************************************************************************
                          Timer.h  -  description
                             -------------------
    begin                : Wed Dec 6 2000
    copyright            : (C) 2000 by Christian Blum
    email                : chr_blum@hotmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef TIMER_H
#define TIMER_H

#include "config.h"

#ifdef _WIN32
#include <winsock.h>
#include <chrono>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include <stdio.h>

class Timer {
private:
#ifdef _WIN32
    chrono::steady_clock::time_point real_time;
#else
  struct rusage res;
  struct timeval tp;
  double virtual_time, real_time;
#endif

public:
  enum TYPE {REAL, VIRTUAL};
  Timer(void);
  double elapsed_time(const TYPE& type);
};
#endif
