/***************************************************************************//**
**  \mainpage Monte Carlo eXtreme - GPU accelerated Monte Carlo Photon Migration \
**      -- OpenCL edition
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2009-2018
**
**  \section sref Reference:
**  \li \c (\b Yu2018) Leiming Yu, Fanny Nina-Paravecino, David Kaeli, and Qianqian Fang,
**         "Scalable and massively parallel Monte Carlo photon transport simulations 
**         for heterogeneous computing platforms," J. Biomed. Optics, 23(1), 010504 (2018)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

#include "tictoc.h"

#define _BSD_SOURCE

#ifndef USE_OS_TIMER

#ifdef MCX_OPENCL

#include <CL/cl.h>
/* use OpenCL timer */
static cl_ulong timerStart, timerStop;
cl_event kernelevent;

unsigned int GetTimeMillis () {
  float elapsedTime;
  clGetEventProfilingInfo(kernelevent, CL_PROFILING_COMMAND_START,
                        sizeof(cl_ulong), &timerStart, NULL);
  clGetEventProfilingInfo(kernelevent, CL_PROFILING_COMMAND_END,
                        sizeof(cl_ulong), &timerStop, NULL);
  elapsedTime=(timerStop - timerStart)*1e-6;
  return (unsigned int)(elapsedTime);
}

unsigned int StartTimer () {
  return 0;
}

#else

#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime_api.h>
/* use CUDA timer */
static cudaEvent_t timerStart, timerStop;

unsigned int GetTimeMillis () {
  float elapsedTime;
  cudaEventRecord(timerStop,0);
  cudaEventSynchronize(timerStop);
  cudaEventElapsedTime(&elapsedTime, timerStart, timerStop);
  return (unsigned int)(elapsedTime);
}

unsigned int StartTimer () {
  cudaEventCreate(&timerStart);
  cudaEventCreate(&timerStop);

  cudaEventRecord(timerStart,0);
  return 0;
}

#endif

#else

static unsigned int timerRes;
#ifndef _WIN32
#if _POSIX_C_SOURCE >= 199309L
#include <time.h>   // for nanosleep
#else
#include <unistd.h> // for usleep
#endif
#include <sys/time.h>
#include <string.h>
void SetupMillisTimer(void) {}
void CleanupMillisTimer(void) {}
long GetTime (void) {
  struct timeval tv;
  timerRes = 1000;
  gettimeofday(&tv,NULL);
  long temp = tv.tv_usec;
  temp+=tv.tv_sec*1000000;
  return temp;
}
unsigned int GetTimeMillis () {
  return (unsigned int)(GetTime ()/1000);
}
unsigned int StartTimer () {
   return GetTimeMillis();
}

#else
#include <windows.h>
#include <stdio.h>
/*
 * GetTime --
 *
 *      Returns the curent time (from some uninteresting origin) in usecs
 *      based on the performance counters.
 */

long GetTime(void)
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

   return ((long) (((double) counter.QuadPart) * cycles_per_usec));
}

#pragma comment(lib,"winmm.lib")

unsigned int GetTimeMillis(void) {
  return (unsigned int)timeGetTime();
}

/*
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
unsigned int StartTimer () {
   SetupMillisTimer();
   return 0;
}
void CleanupMillisTimer(void) {
  if (timeEndPeriod(timerRes) == TIMERR_NOCANDO) {
    fprintf(stderr,"WARNING: bad return value of call to timeEndPeriod.\n");
  }
}

#endif

#endif

#ifdef _WIN32
#include <windows.h>
#elif _POSIX_C_SOURCE >= 199309L
#include <time.h>   // for nanosleep
#else
#include <unistd.h> // for usleep
#endif

/**
  @brief Cross-platform sleep function
*/

void sleep_ms(int milliseconds){
#ifdef _WIN32
    Sleep(milliseconds);
#elif _POSIX_C_SOURCE >= 199309L
    struct timespec ts;
    ts.tv_sec = milliseconds / 1000;
    ts.tv_nsec = (milliseconds % 1000) * 1000000;
    nanosleep(&ts, NULL);
#else
    usleep(milliseconds * 1000);
#endif
}
