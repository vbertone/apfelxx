//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/messages.h"

#include <chrono>

namespace apfel
{
  /**
   * @brief The Timer class computes the time elapsed between start
   * and stop.
   */
  class Timer
  {
  public:
    /**
     * @brief The Timer default constructor.
     */
    Timer() { start(); }

    /**
     * @brief This function starts the timer.
     */
    void start() { startTime = std::chrono::steady_clock::now(); }

    /**
     * @brief This function stops the timer and reports the elapsed
     * time in seconds since the last time the timer was started.
     */
    void stop(bool const& ForceDisplay = false)
    {
      auto end = std::chrono::steady_clock::now();
      auto diff = end - startTime;
      if (GetVerbosityLevel() > 1 || ForceDisplay)
        printf("time elapsed: %5.6f seconds\n", std::chrono::duration <double, std::milli> (diff).count() * 1e-3);
    }

  private:
    std::chrono::time_point<std::chrono::steady_clock> startTime;
  };
}
