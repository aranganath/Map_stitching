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


#include "hough3d.h"
#include "common.h"
#include "io.h"
#include "scan3d.h"

#include <vector>
#include <list>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <fstream>

using namespace h3d;
using namespace std;

const double delta = 0.15;



int main(int argc,char** argv) {

  if ( argc != 3 ) {
    cout << "Please provide the name of two json files"<< endl;
    return 1;
  }

  scan3D scan1;
  scan1.load_scan_jsonfile(argv[1]);

  scan3D scan2;
  scan2.load_scan_jsonfile(argv[2]);  

  scan1.save_scan_to_file("saved_scan");

  /*
  vector<point> points1 = test_shape_pyramid();
  save_vector(points1,"saved_pyramid_points");
  vector<point> points2 = test_shape_pyramid();
  */
  h3d_params p;

  unsigned int i,j;
  double max1 = -1,max2 = -1;
  for ( i = 0 ; i < scan1.get_nscans() ; i++ ) {
    for ( j = 0 ; j < scan1.get_nbeams() ; j++ ) {
      if ( scan1.valid[i][j] )
	max1 = scan1.readings[i][j] > max1 ? scan1.readings[i][j] : max1;
      if ( scan2.valid[i][j] )
	max2 = scan2.readings[i][j] > max2 ? scan2.readings[i][j] : max2;
    }
  }
   
  p.max_norm = 1.1 * ( max1 > max2 ? max1 : max2 );

  h3d_struct ht1(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);
  matrix33 Rt;
  point Tt(0,0,0);

  ht1.compute_cube_ht_normals(scan1,Rt,Tt);
  ht1.save_hs_to_file("saved_spectrum");

  h3d_struct ups = ht1.upsample(p.upsample);
  ups.save_hs_to_file("saved_upsample");

  list<buffer_peak> buffer_peaks = ups.matrix_local_maxima(ups.cube_hs[0]);
  save_list(buffer_peaks,"saved_maxima_face_1");

  buffer_peaks = ups.matrix_local_maxima(ups.cube_hs[1]);
  save_list(buffer_peaks,"saved_maxima_face_2");

  buffer_peaks = ups.matrix_local_maxima(ups.cube_hs[2]);
  save_list(buffer_peaks,"saved_maxima_face_3");

  buffer_peaks = ups.matrix_local_maxima(ups.cube_hs[3]);
  save_list(buffer_peaks,"saved_maxima_face_4");

  buffer_peaks = ups.matrix_local_maxima(ups.cube_hs[4]);
  save_list(buffer_peaks,"saved_maxima_face_5");

  buffer_peaks = ups.matrix_local_maxima(ups.cube_hs[5]);
  save_list(buffer_peaks,"saved_maxima_face_6");
  
  list<s2_peak> sphere_peaks = ups.find_hs_peak();
  save_list(sphere_peaks,"saved_hs_peaks");

  matrix cyl(p.cyl_num_alpha_steps,p.cyl_num_beta_steps);
  cyl = project_hs_on_cylinder(ups,s2_point(0.6973,0.7098,0.996),s2_point(0,0.1390,-0.9903),
			       p.cyl_num_alpha_steps,
			       p.cyl_num_beta_steps,2.6180,true,2);
  
  cyl.save_to_file("saved_cylinder");

  h3d_struct ht2(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);

  matrix33 true_R;
  true_R = rot_axis_angle(s2_point(0.1759,0.6390,-0.7488),60);
  point true_T(0.23,-0.4,0.04);
  
  ht2.compute_cube_ht_normals(scan2,true_R,true_T);

  list<rotation_hypothesis> hyp = ht1.guess_rotation(ht2);
  save_list(hyp,"saved_rotation_hypothesis");
  deque<translation_hypothesis> tran;
  list<roto_translation_hypothesis> roto_tran;
  roto_translation_hypothesis tmp_rototran;

  list<rotation_hypothesis>::iterator r_it;
  for ( r_it = hyp.begin() ; r_it != hyp.end() ; r_it++ ) {
    tran = ht1.guess_translation(ht2,r_it->rotation);
    tmp_rototran.rot = *r_it;
    for ( unsigned int i = 0 ; i < tran.size() ; i++ ) {
      tmp_rototran.tran = tran[i];
      roto_tran.push_back(tmp_rototran);
    }
  }
  
  save_list(roto_tran,"saved_rototranslations");

  vector<ranked_solution> v_final = ht1.evaluate_solutions_hausdorff(scan1,scan2,true_R,true_T,roto_tran);

  save_vector(v_final,"saved_final_results");

  return 0;
	   
}
