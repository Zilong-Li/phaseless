#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>

class Timer
{
protected:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_clock, prev_clock;

public:
    Timer();
    ~Timer();
    void clock();
    unsigned int reltime();
    unsigned int abstime();
};

inline Timer::Timer()
{
    start_clock = std::chrono::high_resolution_clock::now();
}

inline Timer::~Timer()
{
}

inline void Timer::clock()
{
    prev_clock = std::chrono::high_resolution_clock::now();
}

inline unsigned int Timer::reltime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                 prev_clock)
        .count();
}

inline unsigned int Timer::abstime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                 start_clock)
        .count();
}

#endif // TIMER_H_
