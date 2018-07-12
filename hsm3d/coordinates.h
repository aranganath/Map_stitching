 /* AUTORIGHTS
 Copyright (c) 2008 The Regents of the University of California.
 All Rights Reserved.
 
 Created by Stefano Carpin
 University of California, Merced, Robotics Lab - School of Engineering
 
 Permission to use, copy, modify, and distribute this software and its
 documentation for educational, research and non-profit purposes,
 without fee, and without a written agreement is hereby granted,
 provided that the above copyright notice, this paragraph and the
 following three paragraphs appear in all copies.
 
 This software program and documentation are copyrighted by The Regents
 of the University of California. The software program and
 documentation are supplied "as is", without any accompanying services
 from The Regents. The Regents does not warrant that the operation of
 the program will be uninterrupted or error-free. The end-user
 understands that the program was developed for research purposes and
 is advised not to rely exclusively on the program for any reason.
 
 IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
 ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
 CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
 BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
 MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*/

/*!
  \file coordinates.h
  \brief Routines to convert back and forth between the different coordinate systems used.

  Floating point data types are defined though the REAL macro so that at
compile time you can opt for single or double precision. See file precision.h for details.
 */

#ifndef __COORDINATES__
#define __COORDINATES__

//#include <iostream>
#include "common.h"
//#include "precision.h"



// TESTED: OK (see testcoordinates)
/*!
\brief Namespace grouping all functions and classes needed by the scan mathing  algorithms.
 */

namespace h3d {


  /*!
    \brief A point on the face of the cube centered with edge size 2
    and centered in 0,0,0
  */
  struct cube_point {
    /*! Face of the cube. Must be between 1 and 6 */
    unsigned short int face;
    /*! First coordinate on the face. Must be between -1 and 1*/
    REAL u;
    /*! Second coordinate on the face. Must be between -1 and 1*/
    REAL v;
  };

  /*! 
    \brief Buffer representing a discretization of cube_point elements
    as patches on faces
    
    Discretizes the  the cube centered with edge size 2
    and centered in 0,0,0 using the same number of patches on each face.
    Resolution is not stored in the class.
   */
  struct buffer_point {
    /*! Face of the cube. Must be between 1 and 6 */
    unsigned short int face;
    /*! First patch index */
    unsigned int i;
    /*! Second patch index */
    unsigned int j;
  };

  /*!
    \brief A direction in space associated with a scan point.
   */
  struct tt_struct {
    /*! Rotation angle about the y axis*/
    REAL tilt;
    /*! Rotation angle about the z axis*/
    REAL theta;
    /*! Generated from a valid scan point? */
    bool singular;
  };
  
  struct h3d_struct;
  struct matrix;

  /*! 
    \brief Maps a point on the continous cube to a patch on the discretized cube
    @param p point on the continuous cube
    @paarm d resolution used for the discretized cube
   */
  buffer_point coords_cube_to_cell(const cube_point& p,unsigned int d);
  /*!
    \brief Converts a point on the coninuous cube to a poin on the S2 sphere
   */
  s2_point coords_cube_to_s2(const cube_point&);
  /*! \brief Converts a point on S2 to a point on the continuous cube */
  cube_point coords_s2_to_cube(const s2_point&);
  /*! \brief Converts a point on S2 to a point on the discretized cube */
  buffer_point coords_s2_to_cell(const s2_point&,unsigned int);
  /*! \brief Convertes a point on the discretized cube to a point on the continuous cube */
  cube_point coords_cell_to_avg_cube(const buffer_point&,unsigned int);
  void rho2index(const h3d_struct&,REAL,int&,REAL&);
  /*! \brief Projects the Hough spectrum defined on S2 on a surrouding circle */
  matrix project_hs_on_cylinder(const h3d_struct& s,
  		      const s2_point& m,const s2_point& n,
  		      unsigned int alpha_steps,unsigned int beta_steps,REAL amplitude,
  		      bool wba,unsigned int spice);

 
  /*! \brief Convertes the direction associated with a point to the tilt and theta components */
  tt_struct scan_dir2tt(const s2_point&);

  bool operator<(const buffer_point&,const buffer_point&);

}

#endif
