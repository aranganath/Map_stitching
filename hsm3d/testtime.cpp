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


#include "io.h"
#include "hough3d.h"

#include <iostream>
#include <vector>
#include <sys/time.h>
#include <time.h>

using namespace h3d;
using namespace std;

void timer(const char* msg=NULL) {
  static bool initialized = false;
  static struct timeval tv;
  if ( ! initialized ) {
    gettimeofday(&tv,NULL);
    initialized = true;
    return;
  }
  struct timeval tvnow;
  gettimeofday(&tvnow,NULL);

  unsigned int tms = (tvnow.tv_sec * 1000000 + tvnow.tv_usec - tv.tv_sec*1000000 - tv.tv_usec )/ 1000;

  cout << msg << ":" << tms<< " ms" << endl;

  tv = tvnow;
}


int main(int argc,char **argv) {

  if ( argc != 2 ) {
    cerr << "Please provide the filename where points are (tabular form)" << endl;
    return 1;
  }

  h3d_params p;

  vector<point> points1 = load_points_from_table(argv[1]);
  vector<point> points2 = load_points_from_table(argv[1]);
  matrix33 R = rot_z(deg2rad(135)) * rot_y(deg2rad(25));
  apply_rotation(points2,R);
  h3d_struct ht2(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);
  
  ht2.compute_cube_ht(points2);


  cout << "Point clouds consists of " << points1.size() << " points"<<endl;
  
 

  timer();
  h3d_struct ht1(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);
  
  timer("Time to create h3d_struct");

  ht1.compute_cube_ht(points1);

  timer("Time to compute Hough transform");

  h3d_struct ups = ht1.upsample(p.upsample);

  timer("Time to upsample");

  list<s2_peak> sphere_peaks = ups.find_hs_peak();

  timer("Time to find peaks on the Hough spectrum");

  list<rotation_hypothesis> hyp = ht1.guess_rotation(ht2);
  
  timer("Time to guess rotation hypothesis");

  return 0;

}
