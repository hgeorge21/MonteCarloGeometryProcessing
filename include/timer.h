#ifndef TIMER_H
#define TIMER_H

#include <functional>

// Timer for timing a function
// timer function
// copied from https://github.com/libigl/libigl/blob/master/tutorial/717_FastWindingNumber/main.cpp?fbclid=IwAR2sUXck8cvuTrhakA2NGNCvDSXxk04sTMJtUvLtGrZc3vFl_dRJ9TGkn3k
//
// Input:
//   func function to time
// Output:
//   total_time time to run the function

double timer(std::function<void(void)> func);

#endif //TIMER_H
