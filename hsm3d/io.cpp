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

#include <list>
#include <iterator>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>




using namespace std;

namespace h3d {

 // TESTED: WORKS
  ostream& operator<<(ostream& os,const point& p) {

    os << p.x << " " << p.y << " " << p.z;
    return os;

  }

  // TESTED: WORKS
  ostream& operator<<(ostream& os,const point* p) {

    os << p->x << " " << p->y << " " << p->z;
    return os;

  }

 // TESTED: WORKS
   ostream& operator<<(ostream& os,const matrix33& m) {
    
    for (unsigned int i = 0 ; i < 3 ; i++ )
      for (unsigned int j = 0 ; j < 3 ; j++ )
	os << m.data[i][j] << " " ;
    
    return os;
  }

   ostream& operator<<(ostream& os,const cube_point& c) {
    os << c.face << " " << c.u << " " << c.v;
    return os;
  }
  
  ostream& operator<<(ostream& os,const buffer_point& b) {
    os << b.face << " " << b.i << " " << b.j;
    return os;
  }

 // TESTED: works
  ostream& operator<<(ostream& os,const buffer_peak& b) {
    os << b.value << " " << b.p.face << " " << b.p.i << " " << b.p.j;
    return os;
  }
  
  // TESTED: works
  ostream& operator<<(ostream& os,const corr_peak& b) {
    os << b.value << " " << b.pos;
    return os;
  }
  
  // TESTED: works
  ostream& operator<<(ostream& os,const s2_peak& b) {
    os << b.value << " " << b.p.x << " " << b.p.y << " " << b.p.z;
    return os;
  }


  ostream& operator<<(ostream& os,const rotation_hypothesis& h) {

    os << h.value << " " << h.rotation;
    return os;
  }

  ostream& operator<<(ostream& os,const translation_hypothesis& h) {
    
    os << h.value << " " << h.direction;
    return os;
  }


  ostream& operator<<(ostream& os,const roto_translation_hypothesis& h) {
    
    os << h.rot << " " << h.tran;
    return os;
  }

  ostream& operator<<(ostream& os,const ranked_solution& h) {
    
    os << h.value << " " << h.rot << " " << h.tran;
    return os;
  }



  ostream& operator<<(ostream& os,const matrix&M) {

    for ( unsigned int i = 0 ; i < M.r ; i++ )  {
      for ( unsigned int j = 0 ; j < M.c ; j++ )
	os << M.data[i][j] << " ";
      os << endl;
    }
    return os;
  }


  vector<REAL> load_REAL_vector_from_file(const char* fname) {
    deque<REAL> tmp;
    ifstream ifs(fname);
    assert(ifs);
    REAL d;
    while ( ifs ) {
      ifs >> d;
      if ( ! ifs.eof() )
	tmp.push_back(d);
    }
    vector<REAL> retval(tmp.size());
    copy(tmp.begin(),tmp.end(),retval.begin());
    return retval;
  }

  // load points from a file where every row is x,y,z
  vector<point> load_points_from_table(const char* fname) {
    list<point> tmp_list;
    point tmp_point;
    ifstream ifile(fname);
    
    assert(ifile);

    while ( ifile ) {
      ifile >> tmp_point.x;
      ifile >> tmp_point.y;
      ifile >> tmp_point.z;
      tmp_list.push_back(tmp_point);
    }

    vector<point> retval(tmp_list.begin(),tmp_list.end());
    return retval;

  }


  void apply_rotation(vector<point>& v,const matrix33& m) {
    
    vector<point>::iterator it;
    point p;
    for ( it = v.begin() ; it != v.end() ; it++ ) {
      p = *it;
      it->x  = p.x * m.data[0][0] + p.y * m.data[0][1] + p.z *m.data[0][2];
      it->y  = p.x * m.data[1][0] + p.y * m.data[1][1] + p.z *m.data[1][2];
      it->z  = p.x * m.data[2][0] + p.y * m.data[2][1] + p.z *m.data[2][2];
    }
  

  }
  
   
}
