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


#ifndef __HOUGH3D__
#define __HOUGH3D__

//#include "coordinates.h"
//#include "common.h"
#include "scan3d.h"
 #include <string>

//#include "precision.h"

#include <map>
#include <vector>
#include <string>
 // #include <algorithm>
#include <list>
#include <ANN/ANN.h>
//#include <iostream>

/*! \file hough3d.h
  \brief Data types and functions to compute the Hough spectrum in 3 dimensions. 
*/


namespace h3d {

  /*! 
    \brief A peak in a real function befined over the cube of patches
   */
  
  
  struct buffer_peak {
    /*! Location of the peak */
    buffer_point p;
    /*! Value of the peak */
    REAL value;
  };

  /*! \brief A peak in the a real function defined over the unit sphere */
  struct s2_peak {
    /*! Location of the peak*/
    s2_point p;
    /*! Value of the peak */
    REAL value;
  };

  struct obs_obj {

    point M_row;
    REAL Z_row;
    REAL weight;

  };

  /*! \brief A ranked rotation hypotesis */
  struct rotation_hypothesis {
    /*! Rotation hypothesis */
    matrix33 rotation;
    /*! Rotation score */
    REAL value;
  };
  
  /*! \brief A ranked translation hypothesis */
  struct translation_hypothesis {
    /*! Translation/direction hypotesis */
    point direction;
    /*! Translation score */
    REAL value;
  };

  /*! \brief A rototranslation hypothesis.

   This is just a coupling between a rotation and a translation hypothesis. Note that it
does not have its own combined score (yet). */
  struct roto_translation_hypothesis {
    /*! A ranked rotation hypothesis */
    rotation_hypothesis rot;
    /*! A ranked translation hypothesis */
    translation_hypothesis tran;
  };

  /*! \brief A ranked rototranslation hypothesis */
  struct ranked_solution : public roto_translation_hypothesis {
    /*! Actual score */
    REAL value;
    /*! Score of the number of points out of field in the Hausdorff scoring */
    REAL num_out_of_field;
    /*! Score of the number of points hidden in the Hausdorff scoring */
    REAL num_behind_something;
    /*! Score of the number of points not consistent in the Hausdorff scoring */
    REAL num_inconsistent;
    /*! Score of the number of points plausible in the Hausdorff scoring */
    REAL num_plausible;
    
  };
  struct new_ranked_solution : public roto_translation_hypothesis {
    /* The distance between two points is stored (This is analogous to the actual value)*/
    REAL distance;

  };

  /*! \brief View score of a point */
  struct view_score {
    /*! Whether the point is visible */
    bool visible;
    /*! Distance of the point */
    REAL rho;
    /*! Whether the point is at a singular configuration */
    bool singular;
  };
 
 
  bool operator<(const s2_peak&,const s2_peak&);
  bool operator<(const rotation_hypothesis&,const rotation_hypothesis&);
  bool operator<(const ranked_solution&,const ranked_solution&);
  bool operator<(const translation_hypothesis&,const translation_hypothesis&);
  bool operator<(const corr_peak&,const corr_peak&);
  bool operator<(const new_ranked_solution&,const new_ranked_solution&);
  /*! \brief Determines if we already found a rotation hypothesis sufficienlty near to a candidtate one
    @param computed list of already computed rotation hypothesis
    @param newhyp newhypothesis to test
    @param threshold tollerance in degrees
    @return returns true if there exist a rotation hypothesis in the list closer than the given threshold
   */
  bool closer_than(const std::list<rotation_hypothesis>& computed,const rotation_hypothesis& newhyp, int threshold);
  /*! \brief Determines if a given direction is parallel or closer than a certain threshold to a list
    of other directions.
    
   */
  bool closer_than_deg_or_parallel(const std::list<s2_point>& computed,const s2_point& newdir,int threshold);
  
  /*! 
    \brief The actual class storing both the three dimensional Hough transform and Hough
    spectrum.
   */
  class h3d_struct {


  public:
    /*! Set of points during initialization */ 
    std::vector<point> ptsImg;
    /*! Angular resolution in degrees */
    REAL angular_cell_size_deg;
    /*! Distance resolution */
    REAL rho_cell_size;
    /*! Minimum value for rho  */
    REAL rho_min;
    /*! Maximum value for rho */
    REAL rho_max;
    /*! Number of intervals  in which the rho interval is divided */
    unsigned int rho_ncells;
    /*! Number of subdivisions per cube edge (i.e. every face in the cube will have this number squared of patches */
    unsigned int cube_ncells;
    /*! Normals ot the given cube patches. Each cube pacth (i.e. buffe_point) is used to index and retrieve the associated normal */
    std::map<buffer_point,s2_point> cube_normals;
    /*! The actual Hough spectrum in 3D stored in the cube buffer */
    REAL **cube_hs[6];
    /*! The actual Hough transform in 3D stored in the cube buffer */
    REAL ***cube_ht[6];

    /*! An instance of the params class storing the values to use for the solution of the pose registration problem */
    h3d_params params;
  
    /*! \brief Creates an instance with the given resolution
      @param a angular resolution
      @param r rho resolution
      @param rmin minimum value for rho
      @param rmax maximum value for rho
     */
    h3d_struct(REAL a ,REAL r ,REAL rmin,REAL rmax);
    /*! \brief Disposes memory */
    ~h3d_struct();

    /*! \brief Computes the three dimensional Hough transform and spectrum

      This version does not exploits normals (i.e. it is slower) and assumes points are passed in absolute
      coordinates
      @param points vector with absolute coordinates
     */
    void compute_cube_ht(const std::vector<point>& points);

    /*!
      \brief Computes the three dimensional Hough transform and spectrum using the normals

      This version exploits normals and is therefore much faster. Two  additional parameters may be passed to 
      apply a given roto-translation to the three dimensional scan before starting the computation. In an 
      online scenario (i.e. data returned by a robot moving around and repeatedly trying to register successive
      scans), you should simply pass the identity rotation and the null translation.
      @param scan three dimensional scan with the data returned by the scanner
      @param rcorr preliminary rotation correction to apply 
      @param tcorr preliminary translation correction to apply
     */
    void compute_cube_ht_normals(const scan3D& scan,const matrix33& rcorr,const point& tcorr);

    /*! \brief Saves the internal represention to various files. Useful only for debugging */
    void save_hs_to_file(const std::string&);

    /*! \brief Upsamples the transform, thus yielding a smoother function 
      @param order interpolation order
     */
    h3d_struct upsample(unsigned int order) const;
    
    /*! \brief Finds local maxima in a face of the cube */
    std::list<buffer_peak> matrix_local_maxima(REAL **);
    /*! \brief Finds all local maxima in the Hough spectrum (i.e. considers all 6 faces */
    std::list<s2_peak> find_hs_peak();
    /*! \brief Filters away maxima that are to close to each other accordingly to the given resolution. */
    std::list<unsigned int> filter_results(const std::vector<s2_peak>&,REAL);

    /*! \brief Guesses a set of promising rotations to allign two spectra */
    std::list<rotation_hypothesis> guess_rotation(const h3d_struct&);

    /*! \brief Given a candidate rotation hypothesis, guesses a set of promising rotation to allign two spectra */
    std::deque<translation_hypothesis> guess_translation(const h3d_struct&,const matrix33&);

    std::vector<point> choose_sample_points(const scan3D&, unsigned int); 
    
    /*! \brief Evaluates the rototranslation hypothesis and sorts them accordingly  */
    std::vector<ranked_solution> evaluate_solutions_hausdorff(const scan3D&,const scan3D&,const matrix33&,const point&,const std::list<roto_translation_hypothesis>&);
    // std::vector<ranked_solution> evaluate_solutions_hausdorff_abspts(const scan3D&,const scan3D&,const std::list<roto_translation_hypothesis>&);
    std::vector<new_ranked_solution> evaluate_solutions_hausdorff_pts(const std::vector<point>& c, const std::list<roto_translation_hypothesis>& solutions, ANNpointArray &kd_pts, int &npts);
    /*! 
      `\brief Scores a point against a scan
      @param scan 3D scan to score against
      @param o point to score
     */
    view_score scan_whats_your_view(const scan3D& scan,const point& p);
 
    
  private:


    /*! \brief Selects a random subset of points from a scan

      Points are sampled accordingly to their weight, i.e. the area of their associated rectangle, and without reinsertion 
      @param scan 3D scan to sample from
      @param n number of points to sample
      @return vector with sample points. May contain less than n elements
    */
    std::vector<point> sample_points(const scan3D& scan,unsigned int n);

    std::list<std::pair<REAL,REAL> >compensate_one_rotation(const h3d_struct&,const s2_point&,const matrix33&,REAL,unsigned int);


     

  };
  // class h3d_struct_pts{

  //   public:
  //   REAL angular_cell_size_deg;
  //   /*! Distance resolution */
  //   REAL rho_cell_size;
  //   /*! Minimum value for rho  */
  //   REAL rho_min;
  //   /*! Maximum value for rho */
  //   REAL rho_max;




  //   std::vector<point> pointSet;

  //   // Constructor to set all the values for h3d_params
  //   h3d_struct_pts(REAL a ,REAL r ,REAL rmin,REAL rmax);

  //   //Constructor for points in the beginning
  //   h3d_struct_pts(std::vector<point> incoming);

  //   // Constructor to read from a file
  //   h3d_struct_pts(const std::string&);

  //   void compute_cube_ht(const std::vector<point>& c);

  //  };



}


#endif
