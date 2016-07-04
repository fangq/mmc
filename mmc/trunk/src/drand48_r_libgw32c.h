/* Copyright (C) 1995, 1997, 2001 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@gnu.ai.mit.edu>, August 1995.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/***************************************************************************//**
\file    drand48_r_libgw32c.h

\brief   Windows 32 port of drand48_r random number generator from libgw2c
*******************************************************************************/

#ifndef _MMC_POSIX_RAND_LIBGW32C_H
#define _MMC_POSIX_RAND_LIBGW32C_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <ieee754.h>
#include <string.h>

//typedef unsigned long long uint64_t

typedef signed __int8 sint8;
typedef unsigned __int8 uint8;
typedef signed __int16 sint16;
typedef unsigned __int16 uint16;
typedef signed __int32 sint32;
typedef unsigned __int32 uint32;
typedef signed __int64 sint64;
typedef unsigned __int64 uint64;


/** 
\struct drand48_data drand48_r_libgw32c.h
\brief Data structure for communication with thread safe versions.

This type is to be regarded as opaque. It's only exported because users
have to allocate objects of this type.
*/
struct drand48_data
{
  unsigned short int __x[3]; /* Current state. */
  unsigned short int __old_x[3]; /* Old state. */
  unsigned short int __c; /* Additive const. in congruential formula. */
  unsigned short int __init; /* Flag for initializing. */
  uint64 __a; /* Factor in congruential formula. */
};

#ifdef  __cplusplus
extern "C" {
#endif


/* Global state for non-reentrant functions. */
int __drand48_iterate (unsigned short int xsubi[3], struct drand48_data *buffer );
int __erand48_r (unsigned short int xsubi[3],struct drand48_data *buffer, double *result);
int drand48_r (struct drand48_data *buffer, double *result);
int seed48_r (unsigned short int seed16v[3],struct drand48_data *buffer);
double erand48 ( unsigned short int xsubi[3]);


#ifdef  __cplusplus
}
#endif

#endif
