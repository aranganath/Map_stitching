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
  \file common.h
  \brief Math routines and the configuration class.
 */

#ifndef __COMMON__
#define __COMMON__

#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>
#include <vector>
#include <deque>

#include "precision.h"



namespace h3d {

  const REAL EPSILON = 0.00001;

  // forward declarations
  //struct s2_point;

  /*!
    \brief A local maxima in the correlation signal.
   */
  struct corr_peak {
    /*! Maxima value */
    REAL value;
    /*! Maxima position */
    unsigned int pos;
  };

  
  /*!
    \brief A point in three dimensional space.

    Obviously this is used also to represent vectors.
   */
  struct point {
    /* x coordinate */
    REAL x;
    /* y coordinate */
    REAL y;
    /* z coordinate */
    REAL z;

    /*! 
      Creates an instance an initializes all components to the given parameters.
     */
    point(REAL X,REAL Y,REAL Z) : x(X) , y(Y) , z(Z) {};
    /*! 
      Creates an instance an initializes all components to 0.
     */
    point() : x(0),y(0),z(0) {};
   

  };

  /*! Technically an s2_point is a point with length 1.*/
  typedef point s2_point;

  point operator*(REAL,const point&);
  point operator-(const point&,const point&);
  point operator+(const point&,const point&);

  struct matrix33;
  
  /*!
    \brief A matrix with arbitrary dimensions.
   */
  struct matrix {
    /*! The actual matrix data. Public so that you can read write
      without calling setter/getter methods*/
    REAL **data;
    /*! Number of rows */
    unsigned int r;
    /*! Number of columns */
    unsigned int c;
    /*! Creates a matrix of the given size */
    matrix(unsigned int rows,unsigned int columns);
    /*! Creates an instance copying from a matrix33 */
    matrix(const matrix33&);
    /*! Disposes memory */
    ~matrix();
    /*! \brief Normalizes the matrix

      Divides all elements in the matrix by the matrix norm
    */
    void normalize();
    /*! \brief Computes the norm of the matrix.
      
      The norm is defined as the square root of the sum of the squared
      elements of the matrix.
     */
    REAL norm() const;
    /*! Overloaded assignment operator */
    matrix& operator=(const matrix& a);
    /*! Multiplies all elements in the matrix with the given number */
    matrix& operator*=(REAL m);
    /*! Saves the matrix as an ASCII table to a file
     @param filename name of the file */
    void save_to_file(const char* filename) const;
    /*! Resizes the matrix 
      @param newrows new number of rows
      @param newcols new number of columns
     */
    void resize(unsigned int newrows,unsigned int nrecolumns);

  private:
    void allocate_memory();
    void dispose_memory();
    
  };

  /*!
    \brief Multiplies two matrices. 

    Dimensions must be appropriate. Fails (via assert) if not.
   */
  matrix operator*(const matrix& a,const matrix& b);

  /*!
    \brief Subtraction between matrices

    Dimensions must be appropriate. Fails (via assert) if not.
  */
  matrix operator-(const matrix& a,const matrix& b);


  /*! 
    \brief A rotation matrix, i.e. an element of SO(3)
   */
  struct matrix33 {
    /*! 
      \brief Public matirx data data.
       
      This is public to avoid calling setter/getter methods. You should
      never write the elements directly but rather use the appropriate rotation functions.
    */
    REAL data[3][3];
    /*! \brief Initializes the matrix to the identity rotation */
    matrix33();
    /*! \brief  Creates a rotation matrix form a 3x3 generic matrix 

      This method id potentially dangerous because it does not check 
      if the input argument is an element in SO(3).
      @param input matrix to copy from
     */
    matrix33(const matrix& input);
    /*!
      \brief Resets the matrix to the identity transformation.
     */
    void set_identity();
    /*!
      \brief Returns the transpose
     */
    matrix33 transpose() const;
    /*! \brief Overloaded  assignment operator */
    matrix33& operator=(const matrix33&);
  };
  

 
  /*! \brief Multiplies the matrix by a given value. */
  matrix33 operator*(REAL,const matrix33&);
  /*! \brief Returns the product of the given rotation matrix, i.e. the 
    composite rotation. */
  matrix33 operator*(const matrix33&,const matrix33&);
  matrix33 operator+(const matrix33&,const matrix33&);
  /*! 
    \brief Applies a rotation to a given direction vector 
    @param m rotation to apply
    @param d direction to rotate
    @return the rotated direction 
  */
  s2_point operator*(const matrix33& m,const s2_point& d);
  matrix operator*(const matrix&,const matrix&);
  /*!
    \brief Returns a rotation matrix around the x axis
    @param angle rotation around the x axis (in radians)
  */
  matrix33 rot_x(REAL angle);
  /*!
    \brief Returns a rotation matrix around the y axis
    @param angle rotation around the y axis (in radians)
  */
  matrix33 rot_y(REAL agle);
   /*!
    \brief Returns a rotation matrix around the z axis
    @param angle rotation around the z axis (in radians)
  */
  matrix33 rot_z(REAL angle);

  /*!
    \brief Computes the determinant of a rotation matrix
   */
  REAL det33(const matrix33&);
  /*!
    \brief  Computes the inverse of a given rotation matrix
   */
  matrix33 inv33(const matrix33&);

  /*!
    \brief Determines if two given real numbers are closer than  the EPSILON
    constant.
   */
  bool real_compare(const REAL a,const REAL b);

  /*! \brief Converts an angle from degrees to radians */
  REAL deg2rad(REAL angle);
  /*! \brief Converts an angle from radians to degrees */
  int rad2deg(REAL angle);

  // acoordingly to K&R's book,  %'s result is machine dependent  for negative
  // operands, so we need to redefine it properly (sic)
  /*!
    \brief MATLAB compatible reminder between two numbers
    
    According to K&R's book, the result of C's builtin reminder operator
    is machine dependent  for negative
    operands, so it is necessary to properly redefine it.
   */
  int mod(int m,int n); 

  /*!
    \brief Determines the size of the matrix resulting from up-interpolation
    of a given order
    @param r number of rows of the original matrix
    @param c number of columns of the original matrix
    @param order interpolation order
    @param nr number of rows in the interpolated matrix
    @param nc number of columns in the interpolated matrix

  */
  void interp2_size(unsigned int r,unsigned int c,unsigned int order,unsigned int& nr ,unsigned int& nc);

  /*!
    \brief Performs a bidimensional linear interpolation of a given table, idetical to MATLAB's interp2 command
    
    @param src table to interpolate
    @param r number of rows in the starting table
    @param c number of columns in the starting table
    @param order interpolation order
    @param dat resulting table. Space must have been already allocated, accordingly to the results provided by interp2_size
   */
  void interp2(REAL ** src,unsigned int r,unsigned int c,unsigned int order,REAL** dst);

  /*! 
    \brief Determines if the angle between two directions is greater than a given threshold.
    @param a first directions
    @param b second direction
    @param t threshold
  */
  bool apart(const s2_point& a,const s2_point& b,REAL t);

  /*! 
    \brief Computes the cross product between two directions/vectors
   */
  s2_point cross(const s2_point&,const s2_point&);
  /*! \brief Scalar product betwen two directions */
  REAL dot(const s2_point&,const s2_point&);
  /*! \brief Scalar product between two vectors of arbitrary size 
    @param a first vector
    @param b second vectord
    @param n vectors size
   */
  REAL dot_vector(REAL*,REAL*, int);
  /*! \brief Modulus of a point */
  REAL norm(const s2_point&);
  /*! \brief Angle between two directions */
  REAL angle_between_vectors(const s2_point&,const s2_point&);
  /*!
    \bried Builds a rotation matrix from a direction/angle pair.
    @param d direction identifyin the rotation
    @param angle rotation angles
   */
  matrix33 rot_axis_angle(const s2_point& d,REAL angle);
  /*! 
    \brief Returns the direction associated with a vector/point
   */
  s2_point normalize(const point& p);
  /*!
    \brief Computes the circular cross correlation between two vectors of identical size. 

    Important: space for the resulting vector should have been already allocated.
    @param a first vector
    @param b second vector
    @param c correlation (same size as the parameters)
    @param nel number of elements in the correlation vectors
    @param lags vector with indexes of elements in a and b
    @param nlgas number of elements in lags
    
   */
  void circular_cross_correlation(REAL* a,REAL* b,REAL* c, int nel,int* lags, int nlsga);
  /*!
    \brief Computes the cross correlation between two vectors of the same size.
    @param a first vector
    @param b second vector
    @param c resulting vector
    @param lags indices of the result (resized during the computation)
   */
  void cross_correlation(const std::vector<REAL>& a,const std::vector<REAL>& b,std::vector<REAL>& c,std::vector<int>& lags);

  /*!
    \brief Performs a median filtering of a vector, i.e. it replaces each element
    with the median of the n elements surrounding the element itself (n being odd)
    @param a vector to filter
    @param b filgered vector
    @param order number of surrounding elements to consider. At the moment it must be 3 (this will
    be relaxed in a future release.
   */
  void medfilt1(const std::vector<REAL>& a,std::vector<REAL>& out,int order);

  /*!
    \brief Median of three numbers
   */
  REAL med3(REAL a,REAL b,REAL c);

  /*! \brief Pointwise multiplication betwen matrices 
    The two matrices must have the same size. Fails via assert otherwise
  */
  matrix pointwise_prod(const matrix& a,const matrix& b);

  /*!
    \brief Finds all local maxima in a vector of correlation values
    @param v correlation vector
    @param extremes if true the first and the last element may be peaks as well. Deafaults to false
    @return a deque of corr_peak values, i.e. pairs of values and positions.
   */
  std::deque<corr_peak> find_maxima(const std::vector<REAL>& v,bool extremes=false);

  /*!
    \brief Converts a rotation matrix into direction angle representation.
    @param r rotation matrix to convert
    @param d resulting direction
    @param angle resulting angle
   */
  void rot_matrix_to_axis_angle(const matrix33& r,s2_point& d,REAL& angle);

  /*!
    \brief Finds a direction perpendicular to the given one
   */
  s2_point get_one_perpendicular(const s2_point& d);
  
  void solve_least_squares(const std::list<point>&,const std::list<REAL>&,const std::list<REAL>&,matrix&,matrix33&,REAL&);
  
  /*!
    \brief Configuration class for the HSM algorithm.
    
    This is just a collection of paramters.
   */
  struct h3d_params {

    REAL max_norm;
    REAL rho_cell_size;
    REAL angular_cell_size_deg;
    unsigned int upsample;
    unsigned int guess_rot_max_peaks1;
    unsigned int guess_rot_max_peaks2;
    unsigned int guess_rot_first_corr_num_hyp;
    unsigned int guess_rot_max_correction_deg;
    unsigned int guess_rot_num_compensations;
    unsigned int guess_rot_merge_hyp_threshold_deg;
    unsigned int closest_peaks_threshold;
    unsigned int cyl_amplitude_deg;
    unsigned int cyl_num_alpha_steps;
    unsigned int cyl_num_beta_steps;
    bool debug_compensation_show;
    unsigned int guess_tran_max_directions;
    unsigned int guess_tran_max_num_hyp;
    unsigned int guess_tran_two_birds_one_stone;
    unsigned int hausdorff_num_points;
    REAL hausdorff_plausible_threshold;
    bool debug_true_R;
    bool debug_true_T;

    /*! 
      \brief Initializes the parameters to default values.

      Default values are not necessarily optimal. 
     */
    h3d_params();

  };

}

#endif 
