/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2021
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing
**          in Plucker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli,
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a>
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**  \li \c (\b Fang2019) Qianqian Fang and Shijie Yan,
**          <a href="http://dx.doi.org/10.1117/1.JBO.24.11.115002">
**          "Graphics processing unit-accelerated mesh-based Monte Carlo photon transport
**           simulations,"</a> J. of Biomedical Optics, 24(11), 115002 (2019)
**  \li \c (\b Yuan2021) Yaoshen Yuan, Shijie Yan, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-12-1-147">
**          "Light transport modeling in highly complex tissues using the implicit
**           mesh-based Monte Carlo algorithm,"</a> Biomed. Optics Express, 12(1) 147-161 (2021)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

#include "mmc_tictoc.h"

#define _DEFAULT_SOURCE

#ifndef USE_OS_TIMER

#ifdef USE_OPENCL

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
    elapsedTime = (timerStop - timerStart) * 1e-6;
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
    cudaEventRecord(timerStop, 0);
    cudaEventSynchronize(timerStop);
    cudaEventElapsedTime(&elapsedTime, timerStart, timerStop);
    return (unsigned int)(elapsedTime);
}

unsigned int StartTimer () {
    cudaEventCreate(&timerStart);
    cudaEventCreate(&timerStop);

    cudaEventRecord(timerStart, 0);
    return 0;
}

#endif

#else

static unsigned int timerRes;
#ifndef _MSC_VER
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
    gettimeofday(&tv, NULL);
    long temp = tv.tv_usec;
    temp += tv.tv_sec * 1000000;
    return temp;
}
unsigned int GetTimeMillis () {
    return (unsigned int)(GetTime () / 1000);
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

long GetTime(void) {
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
        fprintf(stderr, "Unable to read the performance counter!\n");
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
        fprintf(stderr, "WARNING: Cannot set timer precision.  Not sure what precision we're getting!\n");
    } else {
        timerRes = timeCaps.wPeriodMin;
        fprintf(stderr, "(* Set timer resolution to %d ms. *)\n", timeCaps.wPeriodMin);
    }
}
unsigned int StartTimer () {
    SetupMillisTimer();
    return 0;
}
void CleanupMillisTimer(void) {
    if (timeEndPeriod(timerRes) == TIMERR_NOCANDO) {
        fprintf(stderr, "WARNING: bad return value of call to timeEndPeriod.\n");
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

void sleep_ms(int milliseconds) {
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
