/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/log.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <iomanip> // put_time
#include <sstream>

class Timer
{
  protected:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_clock, prev_clock;

  public:
    Timer();
    ~Timer();
    void clock();
    std::string date();
    double reltime();
    double abstime();
};

inline Timer::Timer()
{
    start_clock = std::chrono::high_resolution_clock::now();
}

inline Timer::~Timer() {}

inline void Timer::clock()
{
    prev_clock = std::chrono::high_resolution_clock::now();
}

inline double Timer::reltime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()
                                                                 - prev_clock)
               .count()
           / 1e3;
}

inline double Timer::abstime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()
                                                                 - start_clock)
               .count()
           / 1e3;
}

inline std::string Timer::date()
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "[%d/%m/%Y-%X]");
    return ss.str();
}

#endif // TIMER_H_
