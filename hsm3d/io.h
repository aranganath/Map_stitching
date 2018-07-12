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


#ifndef __HOUGH_IO
#define __HOUGH_IO

//#include "common.h"
#include "hough3d.h"
//#include "coordinates.h"
//#include "precision.h"

#include <fstream>

//#include <vector>

/*! \file io.h
  \brief  I/O routines to load/save elements from files. Mostly useful for
debudding purposes. 

It also provides overloaded output operators for most data types defined in the library.

*/

namespace h3d{

  /*!
    \brief Loads points for a table in ASCII form.

    Each row must contain space separated x y z values.
    @param filename file containing the values to load
    @return loaded points
   */
  std::vector<point> load_points_from_table(const char* filename);
  void apply_rotation(std::vector<point>&,const matrix33&);

  // g++ 4.2.3 does not recognize this one (but 4.0.3 does :? )
  template<typename T,template<typename> class C>
    void save_container(C<T>& l,const char *fname) {
    std::ofstream of(fname);
    std::ostream_iterator<T> output(of,"\n");
    std::copy(l.begin(),l.end(),output);
  }

  /*! \brief Saves a generic list to a file*/
  template<class T>
    void save_list(std::list<T>& l,const char *fname) {
    std::ofstream of(fname);
    std::ostream_iterator<T> output(of,"\n");
    std::copy(l.begin(),l.end(),output);
  }

  /*! \brief Saves a generic vector to a file */
  template<class T>
    void save_vector(std::vector<T>& l,const char *fname) {
    std::ofstream of(fname);
    std::ostream_iterator<T> output(of,"\n");
    std::copy(l.begin(),l.end(),output);
  }
 
  /*! \brief Saves a generic deque to a file */
  template<class T>
    void save_deque(std::deque<T>& l,const char *fname) {
    std::ofstream of(fname);
    std::ostream_iterator<T> output(of,"\n");
    std::copy(l.begin(),l.end(),output);
  }

  /*! \brief Loads a vector of real numbers from a file. */
  std::vector<REAL> load_REAL_vector_from_file(const char* fn);


  std::ostream& operator<<(std::ostream&,const point&);
  std::ostream& operator<<(std::ostream&,const point*);
  std::ostream& operator<<(std::ostream&,const matrix33&);
  std::ostream& operator<<(std::ostream&,const matrix&);
  std::ostream& operator<<(std::ostream&,const cube_point&);
  std::ostream& operator<<(std::ostream&,const buffer_point&);

  std::ostream& operator<<(std::ostream&,const buffer_peak&);
  std::ostream& operator<<(std::ostream&,const corr_peak&);
  std::ostream& operator<<(std::ostream&,const s2_peak&);
  std::ostream& operator<<(std::ostream&,const rotation_hypothesis&);
  std::ostream& operator<<(std::ostream&,const translation_hypothesis&);
  std::ostream& operator<<(std::ostream&,const roto_translation_hypothesis&);
  std::ostream& operator<<(std::ostream&,const ranked_solution&);
  
}
#endif
