/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
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

/**
   Copyright (C) 1995, 1997, 2001 Free Software Foundation, Inc.
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
   02111-1307 USA.  
*/

/***************************************************************************//**
\file    drand48_r_libgw32c.c

\brief   POSIX 48bit multi-threaded RNG for win32 port via LibGW32C
*******************************************************************************/

#include "drand48_r_libgw32c.h"

struct drand48_data __libc_drand48_data;

int __drand48_iterate (unsigned short int xsubi[3], struct drand48_data *buffer )
{
  uint64 X;
  uint64 result;
  
  /* Initialize buffer, if not yet done. */
  
  // built in bullshit removed
  // if (__builtin_expect (!buffer->__init, 0))
  if (buffer->__init == 0)
  {
  buffer->__a = 0x5deece66dull;
  buffer->__c = 0xb;
  buffer->__init = 1;
  }
  
  /* Do the real work. We choose a data type which contains at least
  48 bits. Because we compute the modulus it does not care how
  many bits really are computed. */
  
  X = (uint64) xsubi[2] << 32 | (uint32) xsubi[1] << 16 | xsubi[0];
  
  result = X * buffer->__a + buffer->__c;
  
  xsubi[0] = result & 0xffff;
  xsubi[1] = (result >> 16) & 0xffff;
  xsubi[2] = (result >> 32) & 0xffff;
  
  return 0;
}

int __erand48_r (unsigned short int xsubi[3],struct drand48_data *buffer, double *result)
{
  union ieee754_double temp;

  /* Compute next state.  */
  if (__drand48_iterate (xsubi, buffer) < 0)
    return -1;

  /* Construct a positive double with the 48 random bits distributed over
     its fractional part so the resulting FP number is [0.0,1.0).  */

  temp.ieee.negative = 0;
  temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
  temp.ieee.mantissa0 = (xsubi[2] << 4) | (xsubi[1] >> 12);
  temp.ieee.mantissa1 = ((xsubi[1] & 0xfff) << 20) | (xsubi[0] << 4);

  /* Please note the lower 4 bits of mantissa1 are always 0.  */
  *result = temp.d - 1.0;

  return 0;
}
//weak_alias (__erand48_r, erand48_r)

int drand48_r (struct drand48_data *buffer, double *result)
{
  return __erand48_r (buffer->__x, buffer, result);
}

int seed48_r (unsigned short int seed16v[3],struct drand48_data *buffer)
{
  /* Save old value at a private place to be used as return value.  */
  memcpy (buffer->__old_x, buffer->__x, sizeof (buffer->__x));

  /* Install new state.  */
  buffer->__x[2] = seed16v[2];
  buffer->__x[1] = seed16v[1];
  buffer->__x[0] = seed16v[0];
  buffer->__a = 0x5deece66dull;
  buffer->__c = 0xb;
  buffer->__init = 1;

  return 0;
}

double erand48 (unsigned short int xsubi[3]){
   double result;
   (void) __erand48_r (xsubi, &__libc_drand48_data, &result);
   return result;
}
   
//weak_alias (__seed48_r, seed48_r)
/*
int main()
{
  struct drand48_data st;
  double t;
  drand48_r(&st,&t);
  printf("%f \n", t );  
  return 0;
}
*/
