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


#include <iostream>
#include <cassert>
#include <cmath>

#include "common.h"
#include "scan3d.h"
#include "hough3d.h"


using namespace std;
using namespace h3d;

int main(int argc,char **argv) {

  if ( argc != 2 ) {
    cout << "Please provide the name of a json fil with a scan"<< endl;
    return 1;
  }

  scan3D scan;
  scan.load_scan_jsonfile(argv[1]);

  h3d_params p;
  h3d_struct ht1(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);

  view_score score;

  double *thetas = new double[scan.get_nbeams()];
  double *tilts = new double[scan.get_nscans()];
  
  unsigned int i,k;
  for ( i = 0 ; i < scan.get_nbeams() ; i++ )
    thetas[i] = deg2rad(scan.get_thetas_deg(i));
  for ( i = 0 ; i < scan.get_nscans() ; i++ ) 
    tilts[i] = deg2rad(scan.get_tilt_deg(i));

  point testpoint = (rot_y(thetas[0]-0.1) * rot_z(tilts[9])) * point(1,0,0);
  
  score = ht1.scan_whats_your_view(scan,testpoint);
  assert( ! score.visible );

  testpoint = (rot_y(thetas[scan.get_nbeams()-1]-0.1) * rot_z(tilts[9])) * point(1,0,0);
  score = ht1.scan_whats_your_view(scan,testpoint);
  assert( ! score.visible );

  testpoint = (rot_y(thetas[9]) * rot_z(tilts[0]-0.1)) * point(1,0,0);
  score = ht1.scan_whats_your_view(scan,testpoint);
  assert( ! score.visible );

  testpoint = (rot_y(thetas[9]) * rot_z(tilts[scan.get_nscans()-1])) * point(1,0,0);
  score = ht1.scan_whats_your_view(scan,testpoint);
  assert( ! score.visible );

  testpoint = point(-1,0,0);
  
  score = ht1.scan_whats_your_view(scan,testpoint);
  assert( ! score.visible );

  
  double theta,tilt;

  for ( k = 0 ; k < scan.get_nscans() ; k++ )
    for ( i = 0 ; i < scan.get_nbeams() ; i++ ) {
      if ( ! scan.valid[k][i] ) {
	cout << "X" ;
	continue;
      }
      theta = deg2rad(scan.get_thetas_deg(i));
      tilt = deg2rad(scan.get_tilt_deg(k));

      testpoint = rot_y(tilt) * rot_z(theta) * point(1,0,0);
      score = ht1.scan_whats_your_view(scan,testpoint);
      assert(score.visible);
      if ( ! score.singular ) {
	assert(fabs(score.rho-scan.readings[k][i])<0.005);
	cout << ".";
      }
      else
	cout << "S";
    }
  cout << endl;

  delete []thetas;
  delete []tilts;

  return 0;

}
