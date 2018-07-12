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
#include "io.h"

#include <list>
#include <iostream>
#include <iterator>
#include <fstream>

using namespace std;
using namespace h3d;

// TESTED: OK

int main(int argc,char **argv) {

  list<s2_point> points;
  //  char filename[255];
  //s2_point p;
  cube_point c_tmp;
  double density = 0.2;
  double u,v;
  unsigned int face;
  /*
  for ( face  = 1 ; face <= 3 ; face++ ) {
    points.clear();
    c_tmp.face = face;
    for ( u = -1 ; u <=1 ; u += density ) {
      c_tmp.u = u;
      for ( v = -1 ; v <=1 ; v += density ) {
	c_tmp.v = v;
	p = coords_cube_to_s2(c_tmp);
	points.push_back(p);
      }
    }
    sprintf(filename,"points_face_%d",face);
    save_container(points,filename);

  }
  */

  cout << "Testing coordinate changes from cube to s2...." << endl;

  s2_point s;
  cube_point c;

  for ( face = 1 ; face <= 6 ; face++ ) {
    for ( u = -1 + density ; u <= 1 - density ; u+= 0.1 ) {
      for ( v = -1 + density ; v <= 1 - density ; v+= 0.1 ) {
	c_tmp.face = face;
	c_tmp.u = u;
	c_tmp.v = v;
	s = coords_cube_to_s2(c_tmp);
	c = coords_s2_to_cube(s);
	if ( ( c.face != c_tmp.face )  ||
	     ( ! h3d::real_compare(c.u,c_tmp.u) ) ||
	     ( ! h3d::real_compare(c.v,c_tmp.v) ) )
	  cout << "Anomaly detected: " <<  c_tmp << " " << c << endl;
      }
    }
  }
  cout << "Done" <<endl;
  cout << "Testing coordinate changes from cube to buffer...." << endl;
  unsigned int width = 10;
  buffer_point bp,bp2;
  for ( bp.face = 1 ; bp.face <= 6 ; bp.face++ ) 
    for ( bp.i = 1 ; bp.i <= width ; bp.i++ )
      for ( bp.j = 1 ; bp.j <= width ; bp.j++ ) {
	c = coords_cell_to_avg_cube(bp,width);
	bp2 = coords_cube_to_cell(c,width);
	if ( ( bp.face != bp2.face )  ||
	     ( bp.i != bp2.i ) ||
	     ( bp.j != bp2.j ) )
	  cout << "Anomaly detected: " <<  bp << " " << bp2 << endl;
      }
  cout << "Done" <<endl;
  return 0;

}
