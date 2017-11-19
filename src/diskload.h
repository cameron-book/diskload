/****************************************************************
 * Copyright (C) 2017 by Michael Bevis, Jim Fowler, Crichton Ogle
 *
 * This file is part of diskload
 *
 * diskload is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *                                                                          
 * diskload is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                          
 * You should have received a copy of the GNU Lesser General Public
 * License.  If not, see http://www.gnu.org/licenses/
 *
 */

/**
 * @file diskload.h
 * @author Michael Bevis, Jim Fowler, Crichton Ogle
 * @date 18 Nov 2017
 * @brief Various methods for computing the elastic response to a disk load
 *
 * @see https://academic.oup.com/gji/article-abstract/205/3/1804/657155
 */

#ifndef _DISKLOAD_H
#define _DISKLOAD_H

/*----------------------------------------------------------------*/

/**
 * @brief Switch for two versions of spherical disk load computations
 *
 * The disk load can be defined in two ways, differing in how mass
 * conservation is addressed.  In the uncompensated case, the loading
 * function for a disk of angular radius \f$\alpha\f$ is
 * \f[
 *   \sigma(\vartheta) = \rho \begin{cases} T, & \mbox{if } 0 \leq \vartheta \leq \alpha \\ 0, & \mbox{if not.} \end{cases}
 * \f]
 * In the compensated case, the loading function is instead
 * \f[
 *   \sigma(\vartheta) = \rho \begin{cases} T, & \mbox{if } 0 \leq \vartheta \leq \alpha \\ -\left( \frac{1 - \cos \alpha}{1 + \cos \alpha} \right) T, & \mbox{if not.} \end{cases}
 * \f] 
 */
typedef enum {
  Uncompensated = 0, /**< The load function is constant inside the disk, and zero outside the disk. */
  Compensated = 1 /**< The load function is constant inside the disk and outside the disk, but arranged outside to ensure mass conservation */
} DiskLoadType;

/*----------------------------------------------------------------*/

/**
 * @brief spherical Earth model parameters
 */
typedef struct {
  double radius;   /**< radius of the Earth in km. */
  double gravity; /**< surface gravity in m/s^2. */
} EarthModel;

/**
 * @brief  Default values for the radius of earth and its surface gravity
 * @relates EarthModel
 */
extern const EarthModel DefaultEarthModel;

/*----------------------------------------------------------------*/

typedef enum 
{
    E_SUCCESS = 0,
    E_LOVE_NUMBER_VECTOR_TOO_SHORT = -1,
    E_LOVE_NUMBER_L0_AND_K0_ARE_NONZERO = -2
} DiskLoadError;

/*----------------------------------------------------------------*/

/**
 * @brief elastic loading Love numbers
 *
 * Most users will allocate `LoveNumbers` with a call to read_love_numbers().
 *
 */
typedef struct {
  int degrees; /**< the length of this vector */
  double* h;  /**< `h[i]` is the Love number \f$h_i\f$ with `i` between 0 and `degrees-1` */
  double* l;  /**< `l[i]` is the Love number \f$\ell_i\f$ with `i` between 0 and `degrees-1` */
  double* k;  /**< `k[i]` is the Love number \f$k_i\f$ with `i` between 0 and `degrees-1` */
} LoveNumbers;

LoveNumbers* diskload_read_love_numbers(const char* filename);

DiskLoadError diskload_truncated(double alpha,DiskLoadType icomp,double theta,double w,int nmax,LoveNumbers* LN,const EarthModel* EM, double *u, double *v, double *g);

void diskload_extrapolate_love_numbers( LoveNumbers* love, int nmax );

void diskload_perror( const char* str );

#endif /* _DISKLOAD_H */
