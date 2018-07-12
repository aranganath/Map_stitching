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


#include "coordinates.h"
#include "common.h"
#include "hough3d.h"
#include "io.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>

using namespace std;

// Functions tested with testcoordinates. Appear to work as of 8/12/2008

namespace h3d {

  // TESTED: WORKS
  buffer_point coords_cube_to_cell(const cube_point& src,unsigned int length) {

    assert(length > 0);
    assert((src.face >=1 ) && (src.face <=6 ));
    assert((src.u >= -1 ) && ( src.u <= 1));
    assert((src.v >= -1 ) && ( src.v <= 1));

    buffer_point retval;

    REAL u = (src.u + 1)/2;
    REAL v = (src.v + 1)/2;

    unsigned int n = length;

    retval.face = src.face;

    retval.i = real_compare(u,1) ? n : 1+floor(u*n);
    retval.j = real_compare(v,1) ? n : 1+floor(v*n);

    return retval;
    
  }

  // TESTED: WORKS
  s2_point coords_cube_to_s2(const cube_point& src) {

    assert((src.face >=1 ) && (src.face <=6 ));
    assert((src.u >= -1 ) && ( src.u <= 1));
    assert((src.v >= -1 ) && ( src.v <= 1));    
    
    s2_point retval;

    switch ( src.face ) {
    case 1:
      retval.x = 1;
      retval.y = src.u;
      retval.z = src.v;
      break;
    case 4:
      retval.x = -1;
      retval.y = -src.u;
      retval.z = -src.v;
      break;
    case 2:
      retval.x = src.u;
      retval.y = 1;
      retval.z = src.v;
      break;
    case 5:
      retval.x = -src.u;
      retval.y = -1;
      retval.z = -src.v;
      break;
    case 3:
      retval.x = src.u;
      retval.y = src.v;
      retval.z = 1;
      break;
    case 6:
      retval.x = -src.u;
      retval.y = -src.v;
      retval.z = -1;
      break;
    }

    REAL norm = sqrt(retval.x*retval.x + retval.y*retval.y + retval.z*retval.z);
    retval.x /= norm;
    retval.y /= norm;
    retval.z /= norm;

    return retval;
  }

  // TESTED: WORKS
  cube_point coords_s2_to_cube(const s2_point& s1) {

    cube_point retval;
    s2_point s = s1;
   
    REAL s_norm =sqrt(s.x*s.x + s.y*s.y + s.z*s.z);
    s.x = s.x / s_norm;
    s.y = s.y / s_norm;
    s.z = s.z / s_norm;

    s2_point sa(fabs(s.x), fabs(s.y), fabs(s.z));

    if ( ( sa.x >= sa.y ) && ( sa.x >= sa.z ) ) {
      retval.face = s.x > 0 ? 1 : 4;
      retval.u = s.y / s.x;
      retval.v = s.z / s.x;
      return retval;
    }
    
    if ( ( sa.y >= sa.x ) && ( sa.y >= sa.z ) ) {
      retval.face = s.y > 0 ? 2 : 5;
      retval.u = s.x / s.y;
      retval.v = s.z / s.y;
      return retval;
    }

    if ( ( sa.z >= sa.y ) && ( sa.z >= sa.x ) ) {
      retval.face = s.z > 0 ? 3 : 6;
      retval.u = s.x / s.z;
      retval.v = s.y / s.z;
      return retval;
    }

    cerr << "Error in coords_s2_to_cube"<< endl;
    cout << "Error point " << s1 << endl;
    exit(1);
    
  }

  // TESTED: WORKS
  buffer_point coords_s2_to_cell(const s2_point& s,unsigned int width) {
    cube_point tmp = coords_s2_to_cube(s);
    buffer_point retval = coords_cube_to_cell(tmp,width);
    return retval;
  }

  // TESTED: WORKS
  cube_point coords_cell_to_avg_cube(const buffer_point& src,unsigned int width) {
     assert((src.face >=1 ) && (src.face <=6 ));
     assert((src.i >= 1 ) && ( src.i <= width));
     assert((src.j >= 1 ) && ( src.j <= width));  

     cube_point retval;
     retval.face = src.face;
     retval.u = (src.i - 0.5 ) / (static_cast<REAL>(width));
     retval.u = retval.u*2 - 1;

     retval.v = (src.j - 0.5 ) / (static_cast<REAL>(width));
     retval.v = retval.v*2 - 1;

     return retval;   
  }

  // TESTED: WORKS
  void rho2index(const h3d_struct& s,REAL rho,int& r1,REAL& alpha) {

    if ( ( rho < s.rho_min ) || ( rho > s.rho_max ) ) {
      r1 = 0;
      alpha = 0;
    }
    else {
      REAL x = s.rho_ncells * ( ( rho - s.rho_min) / (s.rho_max - s.rho_min) );
      r1 = static_cast<int>(floor(x));
      alpha = ( r1 + 0.5 ) - x;

      // useless cast just to get rid of a warning
      assert( (r1>=0 ) && ( r1 < static_cast<int>(s.rho_ncells) ) );
    }

  }

  // Important: memory for res should have been already allocated
  // TESTED: WORKS
  matrix project_hs_on_cylinder(const h3d_struct&s,
			      const s2_point& m,const s2_point& n,
			      unsigned int num_alpha_steps,
			      unsigned int num_beta_steps,
			      REAL amplitude,
			      bool weight_by_alpha,
			      unsigned int spice) {

    matrix res(num_alpha_steps,num_beta_steps);

    unsigned int a,b;
    REAL alpha,sin_alpha,beta;
    s2_point direction;
    buffer_point bp;
    for ( a = 0 ; a < num_alpha_steps ; a++ )
      fill(res.data[0],res.data[0]+num_beta_steps,0);

    for ( a = 0 ; a < num_alpha_steps ; a++ ) {

      alpha = ( M_PI-amplitude)/2 + amplitude * a / (num_alpha_steps -1 ) + deg2rad(spice*2);
      sin_alpha = weight_by_alpha ? sinf(alpha) : 1;
      for ( b = 0 ; b < num_beta_steps ; b++ ) {
	beta = 2 * M_PI * b / (num_beta_steps-1);
	direction = ( rot_axis_angle(m,beta) * rot_axis_angle(n,alpha) ) * m;
	bp = coords_s2_to_cell(direction,s.cube_ncells);
	res.data[a][b] = s.cube_hs[bp.face-1][bp.i-1][bp.j-1] * sin_alpha;
      }
    }

    return res;
    
  }

  // TESTED: WORKS
  tt_struct scan_dir2tt(const s2_point& p) {
    tt_struct r;

    r.singular = false;
    
    if ( ( fabs(p.z) < 0.000001 ) && ( fabs(p.x) <  0.000001 ) ) {
      r.singular = true;
      r.tilt = 0;
    }
    else 
      r.tilt = -atan2(p.z,p.x);

    s2_point p2d = rot_y(-r.tilt)*p;
    assert(fabs(p2d.z) <  0.0001 );
    r.theta = atan2(p2d.y,p2d.x);

    return r;

  }

 
  // TESTED: WORKS
  bool operator<(const buffer_point& a,const buffer_point& b) {

    if ( a.face != b.face )
      return (a.face < b.face);
    if ( a.i != b.i )
      return (a.i < b.i );
    return a.j < b.j;

  }

}
