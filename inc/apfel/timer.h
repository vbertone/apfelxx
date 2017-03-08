//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

namespace apfel
{
  /**
   * @brief The Timer class.
   *
   * Computes the calculation time.
   */
  class Timer {

   public:
    //! Starts the timer.
    void start(){ gettimeofday(&startTime, NULL); }

    //! Stops the timer.
    double stop()
    {
      timeval endTime;
      long seconds, useconds;
      double duration;
      gettimeofday(&endTime, NULL);
      seconds  = endTime.tv_sec  - startTime.tv_sec;
      useconds = endTime.tv_usec - startTime.tv_usec;
      duration = seconds + useconds/1E6f;
      return duration;
    }

    /*!
     * \brief Prints enlapsed time
     * \param duration input time from Timer::stop()
     */
    static void printTime(double const& duration)
    {
      printf("elapsed time: %5.6f seconds\n", duration);
    }

  private:
   timeval startTime;
  };

}
