/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    tictoc.c

\brief   Timing functions for different platforms with libc, CUDA and OpenCL
*******************************************************************************/

#include "tictoc.h"

#ifndef USE_OS_TIMER

#ifdef MCX_OPENCL

#include <CL/cl.h>
/* use OpenCL timer */
static cl_ulong timerStart, timerStop;
cl_event kernelevent;

/**
 * @brief OpenCL timing function using clGetEventProfilingInfo
 *
 * Use OpenCL events to query elapsed time in ms between two events
 */

unsigned int GetTimeMillis () {
  float elapsedTime;
  clGetEventProfilingInfo(kernelevent, CL_PROFILING_COMMAND_START,
                        sizeof(cl_ulong), &timerStart, NULL);
  clGetEventProfilingInfo(kernelevent, CL_PROFILING_COMMAND_END,
                        sizeof(cl_ulong), &timerStop, NULL);
  elapsedTime=(timerStop - timerStart)*1e-6;
  return (unsigned int)(elapsedTime);
}

/**
 * @brief Start OpenCL timer
 */

unsigned int StartTimer () {
  return 0;
}

#else                          /**< use CUDA time functions */

#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime_api.h>
/* use CUDA timer */
static cudaEvent_t timerStart, timerStop;

/**
 * @brief CUDA timing function using cudaEventElapsedTime
 *
 * Use CUDA events to query elapsed time in ms between two events
 */

unsigned int GetTimeMillis () {
  float elapsedTime;
  cudaEventRecord(timerStop,0);
  cudaEventSynchronize(timerStop);
  cudaEventElapsedTime(&elapsedTime, timerStart, timerStop);
  return (unsigned int)(elapsedTime);
}

/**
 * @brief Start CUDA timer
 */

unsigned int StartTimer () {
  cudaEventCreate(&timerStart);
  cudaEventCreate(&timerStop);

  cudaEventRecord(timerStart,0);
  return 0;
}

#endif

#else                          /**< use host OS time functions */

static unsigned int timerRes;
#ifndef MSVC
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
void SetupMillisTimer(void) {}
void CleanupMillisTimer(void) {}

/**
 * @brief Unix timing function using gettimeofday
 *
 * Use gettimeofday to query elapsed time in us between two events
 */

unsigned long GetTime (void) {
  struct timeval tv;
  timerRes = 1000;
  gettimeofday(&tv,NULL);
  unsigned long temp = tv.tv_usec;
  temp+=tv.tv_sec*1000000;
  return temp;
}

/**
 * @brief Convert us timer output to ms
 */

unsigned int GetTimeMillis () {
  return (unsigned int)(GetTime ()/1000);
}

/**
 * @brief Start UNIX timer
 */

unsigned int StartTimer () {
   return GetTimeMillis();
}

#else
#include <windows.h>
#include <stdio.h>

/**
 * @brief Windows timing function using QueryPerformanceCounter
 *
 * Retrieves the current value of the performance counter, which is a high 
 * resolution (<1us) time stamp that can be used for time-interval measurements.
 */

unsigned long GetTime(void)
{
   static double cycles_per_usec;
   LARGE_INTEGER counter;

   if (cycles_per_usec == 0) {
      static LARGE_INTEGER lFreq;
      if (!QueryPerformanceFrequency(&lFreq)) {
         fprintf(stderr, "Unable to read the performance counter frquency!\n");
         return 0;
      }

      cycles_per_usec = 1000000 / ((double) lFreq.QuadPart);
   }

   if (!QueryPerformanceCounter(&counter)) {
      fprintf(stderr,"Unable to read the performance counter!\n");
      return 0;
   }

   return ((unsigned long) (((double) counter.QuadPart) * cycles_per_usec));
}

#pragma comment(lib,"winmm.lib")

unsigned int GetTimeMillis(void) {
  return (unsigned int)timeGetTime();
}

/**
  @brief Set timer resolution to milliseconds
  By default in 2000/XP, the timeGetTime call is set to some resolution
  between 10-15 ms query for the range of value periods and then set timer
  to the lowest possible.  Note: MUST make call to corresponding
  CleanupMillisTimer
*/

void SetupMillisTimer(void) {

  TIMECAPS timeCaps;
  timeGetDevCaps(&timeCaps, sizeof(TIMECAPS)); 

  if (timeBeginPeriod(timeCaps.wPeriodMin) == TIMERR_NOCANDO) {
    fprintf(stderr,"WARNING: Cannot set timer precision.  Not sure what precision we're getting!\n");
  }else {
    timerRes = timeCaps.wPeriodMin;
    fprintf(stderr,"(* Set timer resolution to %d ms. *)\n",timeCaps.wPeriodMin);
  }
}

/**
  @brief Start system timer
*/

unsigned int StartTimer () {
   SetupMillisTimer();
   return 0;
}

/**
  @brief Reset system timer
*/

void CleanupMillisTimer(void) {
  if (timeEndPeriod(timerRes) == TIMERR_NOCANDO) {
    fprintf(stderr,"WARNING: bad return value of call to timeEndPeriod.\n");
  }
}

#endif

#endif
