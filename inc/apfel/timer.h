//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/messages.h"

#include <iostream>
#include <chrono>

using namespace std;

namespace apfel
{
  /**
   * @brief The Timer class.
   *
   * Computes the calculation time.
   */
  class Timer
  {
  public:
    // Constructor
    Timer() { start(); }

    // Starts the timer.
    void start() { startTime = chrono::steady_clock::now(); }

    // Stops the timer.
    void stop()
    {
      auto end = chrono::steady_clock::now();
      auto diff = end - startTime;
      if (GetVerbosityLevel() > 1)
	printf("time elapsed: %5.6f seconds\n", chrono::duration <double, milli> (diff).count() * 1e-3);
    }

  private:
    chrono::time_point<chrono::steady_clock> startTime;
  };
}
