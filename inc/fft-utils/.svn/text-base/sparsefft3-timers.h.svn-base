/* Copyright (C) 2000 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/********************* Timer code (from FFTW) **************************/

#define SPARSE3DFFT_ENABLE_PENTIUM_TIMER 0
#define SPARSE3DFFT_CYCLES_PER_SEC 550000000

#define HAVE_GETTIMEOFDAY 1

#if defined(__GNUC__) && defined(__i386__) && SPARSE3DFFT_ENABLE_PENTIUM_TIMER

/*
 * Use internal Pentium register (time stamp counter). Resolution
 * is 1/SPARSEFFT3_CYCLES_PER_SEC seconds (e.g. 5 ns for Pentium 200 MHz).
 * (This code was contributed by Wolfgang Reimer)
 */

typedef unsigned long long sparse3DFFT_time;

static __inline__ sparse3DFFT_time read_tsc( )
{
  sparse3DFFT_time ret;

  __asm__ __volatile__("rdtsc": "=A" (ret));
  /* no input, nothing else clobbered */
  return ret;
}

#define sparse3DFFT_get_time( )  read_tsc( )
#define sparse3DFFT_time_diff( t1, t2 ) ( ( t1 ) - ( t2 ) )
#define sparse3DFFT_time_to_sec( t ) ( ( ( double ) ( t ) ) / SPARSE3DFFT_CYCLES_PER_SEC )

#define SPARSE3DFFT_TIME_MIN ( 1.0e-4 )

#elif defined(CRAY) || defined(_UNICOS) || defined(_CRAYMPP)

extern double SECONDR( void );

typedef double sparse3DFFT_time;

#define sparse3DFFT_get_time( ) SECONDR( )
#define sparse3DFFT_time_diff( t1, t2 ) ( ( t1 ) - ( t2 ) )
#define sparse3DFFT_time_to_sec( t ) ( t )

#define SPARSE3DFFT_TIME_MIN ( 1.0e-1 )

#elif HAVE_GETTIMEOFDAY

#include <sys/time.h>
#include <unistd.h>

typedef struct timeval sparse3DFFT_time;

extern sparse3DFFT_time sparse3DFFT_gettimeofday_get_time( void );
extern sparse3DFFT_time sparse3DFFT_gettimeofday_time_diff( sparse3DFFT_time t1, sparse3DFFT_time t2 );

#define sparse3DFFT_get_time( ) sparse3DFFT_gettimeofday_get_time( )
#define sparse3DFFT_time_diff( t1, t2 ) sparse3DFFT_gettimeofday_time_diff( t1, t2 )
#define sparse3DFFT_time_to_sec( t ) ( ( double )( t ).tv_sec + ( double )( t ).tv_usec * 1.0E-6 )

#define SPARSE3DFFT_TIME_MIN ( 1.0e-2 )

#define USING_GETTIMEOFDAY 1

extern sparse3DFFT_time sparse3DFFT_gettimeofday_get_time( void );
extern sparse3DFFT_time sparse3DFFT_gettimeofday_time_diff( sparse3DFFT_time t1, sparse3DFFT_time t2 );

#else /* default to clock() */

#include <time.h>

typedef clock_t sparse3DFFT_time;

#ifndef CLOCKS_PER_SEC
#ifdef sun
/* stupid sunos4 prototypes */
#define CLOCKS_PER_SEC 1000000
extern long clock(void);
#else                           /* not sun, we don't know CLOCKS_PER_SEC */
#error Please define CLOCKS_PER_SEC
#endif
#endif

#define sparse3DFFT_get_time( ) clock( )
#define sparse3DFFT_time_diff( t1, t2 ) ( ( t1 ) - ( t2 ) )
#define sparse3DFFT_time_to_sec( t ) ( ( ( double ) ( t ) ) / CLOCKS_PER_SEC )

#define SPARSE3DFFT_TIME_MIN ( 2.0e-1 )

#endif

#define SPARSE3DFFT_TIME_REPEAT 4
#define SPARSE3DFFT_TIME_LIMIT 2.0

#define TIME_FFT( fft, t, tmin ) \
{ \
     sparse3DFFT_time ts,te; \
     double real_total_t = 0, total_t, tcur, tbest = 1e20; \
     int tfft_iters = 1, tfft_iter, tfft_repeat; \
     for ( tfft_repeat = 0; tfft_repeat < SPARSE3DFFT_TIME_REPEAT; ++tfft_repeat ) { \
     do { \
          ts = sparse3DFFT_get_time( ); \
          for ( tfft_iter = 0; tfft_iter < tfft_iters; ++tfft_iter ) { fft; } \
          te = sparse3DFFT_get_time( ); \
          tcur = ( total_t = sparse3DFFT_time_to_sec( sparse3DFFT_time_diff( te, ts ) ) ) / tfft_iters; \
          tfft_iters *= 2; \
     } while ( total_t < tmin ); \
     real_total_t += total_t; \
     if ( tcur < tbest ) tbest = tcur; \
     if ( real_total_t > SPARSE3DFFT_TIME_LIMIT ) break; \
     } \
     t = tbest; \
}

