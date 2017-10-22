//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <iostream>
#include <chrono>

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
    //! Starts the timer.
    void start(){  startTime = std::chrono::steady_clock::now(); }

    //! Stops the timer.
    void stop()
    {
      auto end = std::chrono::steady_clock::now();
      auto diff = end - startTime;

      printf("time elapsed: %5.6f seconds\n",
             std::chrono::duration <double, std::milli> (diff).count()*1e-3);
    }

  private:
    std::chrono::time_point<std::chrono::steady_clock> startTime;
  };
}
