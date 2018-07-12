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


#include "common.h"
#include "io.h"

#include <cmath>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cassert>



using namespace std;

namespace h3d{

  
  // TESTED: WORKS
  bool real_compare(const REAL a, const REAL b) {

    return (fabs(a-b)<EPSILON);
  
  }
  
  // TESTED: WORKS
  REAL deg2rad(REAL deg) {

    return ( deg * 2*M_PI / 360 );

  }

  // TESTED: WORKS
  int rad2deg(REAL rad) {

    return round(360*rad/(2*M_PI));

  }

  // compatible with MATLAB's definition
  int mod(int x,int y) {
    if ( y == 0 ) return x;
    if ( x == y ) return 0;
    int n = floor(static_cast<REAL>(x)/y);
    int r = x - n*y;
    return ( r * y ) < 0 ? -r : r;
  }

 

  // TESTED: WORKS   
  bool apart(const s2_point& a,const s2_point& b,REAL t) {
    
    REAL dp = dot(a,b);
    if ( dp > 1 ) // get rid of rouding errors...
      dp = 1;
    else if ( dp < -1 )
      dp = -1;
    REAL angle =  acos(dp);
    return ( fabs(angle) > deg2rad(t) ) ? true : false;

  }



  point operator-(const point&a,const point&b) {
    point r;
    r.x = a.x-b.x;
    r.y = a.y-b.y;
    r.z = a.z-b.z;
    return r;
  }

  point operator+(const point&a,const point&b) {
    point r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    return r;
  }

  point operator*(REAL f,const point& p) {
    point r(f*p.x,f*p.y,f*p.z);
    return r;
  }

  // TESTED: WORKS   
  s2_point cross(const s2_point& a,const s2_point& b) {

    s2_point r;

    r.x = a.y*b.z - a.z*b.y;
    r.y = -a.x*b.z + a.z*b.x;
    r.z = a.x*b.y - a.y*b.x;

    return r;

  }


  // TESTED: WORKS   
  REAL dot(const s2_point&a,const s2_point&b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
  }

  REAL norm(const s2_point& a) {
    return sqrt(dot(a,a));
  }

  // TESTED: WORKS   
  REAL angle_between_vectors(const s2_point&a,const s2_point&b) {
    return acosf(dot(a,b)/(norm(a)*norm(b)));
  }


  matrix33::matrix33() {
    for (unsigned int i = 0 ; i < 3 ; i++ ) 
      for ( unsigned int j = 0 ; j < 3 ; j++ )
	data[i][j]=0;

    data[0][0] = data[1][1] = data[2][2] = 1;
  }

  matrix33::matrix33(const matrix& m) {
    assert((m.r == 3) && ( m.c == 3 ) );
    for (unsigned int i = 0 ; i < 3 ; i++ ) 
      for ( unsigned int j = 0 ; j < 3 ; j++ )
	data[i][j]=m.data[i][j];
    
  }

  void matrix33::set_identity() {
    for (unsigned int i = 0 ; i < 3 ; i++ ) 
      for ( unsigned int j = 0 ; j < 3 ; j++ )
	data[i][j]=0;
    data[0][0] = data[1][1] = data[2][2] = 1;
  }

  matrix33 operator*(REAL t,const matrix33& m) {
    matrix33 r;

    for (unsigned int i = 0 ; i < 3 ; i++ ) 
      for ( unsigned int j = 0 ; j < 3 ; j++ )
	r.data[i][j] = m.data[i][j] * t;

    return r;
  }

  // returns the determinant of a 3x3 matrix
  // TESTED: WORKS
  REAL det33(const matrix33& m) {
    REAL r;
    r = m.data[0][0]*(m.data[1][1]*m.data[2][2]-m.data[1][2]*m.data[2][1]);
    r-= m.data[0][1]*(m.data[1][0]*m.data[2][2]-m.data[1][2]*m.data[2][0]);
    r+= m.data[0][2]*(m.data[1][0]*m.data[2][1]-m.data[1][1]*m.data[2][0]);
    return r;
  }

  // returns the inverse of a 3x3 matrix. No iterations involved
  // TESTED: WORKS
  matrix33 inv33(const matrix33&m) {

    
    matrix33 r;
    REAL det = 1/det33(m);
    assert(det != 0 );
    
    r.data[0][0] = m.data[1][1]*m.data[2][2]-m.data[1][2]*m.data[2][1];
    r.data[1][0] = -1*(m.data[1][0]*m.data[2][2]-m.data[1][2]*m.data[2][0]);
    r.data[2][0] = m.data[1][0]*m.data[2][1]-m.data[1][1]*m.data[2][0];
    
    r.data[0][1] = -1*(m.data[0][1]*m.data[2][2]-m.data[0][2]*m.data[2][1]);
    r.data[1][1] = m.data[0][0]*m.data[2][2]-m.data[0][2]*m.data[2][0];
    r.data[2][1] = -1*(m.data[0][0]*m.data[2][1]-m.data[0][1]*m.data[2][0]);
    
    r.data[0][2] = m.data[0][1]*m.data[1][2]-m.data[0][2]*m.data[1][1];
    r.data[1][2] = -1*(m.data[0][0]*m.data[1][2]-m.data[0][2]*m.data[1][0]);
    r.data[2][2] = m.data[0][0]*m.data[1][1]-m.data[0][1]*m.data[1][0];

    r = det*r;

    return r;

  }

  // TESTED: WORKS
  matrix33 operator*(const matrix33&a,const matrix33&b) {
    matrix33 r;
    unsigned int i,j,k;

    for ( i = 0 ; i < 3 ; i++ ) 
      for (  j = 0 ; j < 3 ; j++ ) {
	r.data[i][j] = 0;
	for ( k = 0 ; k < 3 ; k++ )
	  r.data[i][j] += a.data[i][k]*b.data[k][j];
      }

    return r;
  }

  matrix33 matrix33::transpose() const {
    matrix33 retval = *this;

    retval.data[0][1] = data[1][0];
    retval.data[0][2] = data[2][0];
    retval.data[1][0] = data[0][1];
    retval.data[1][2] = data[2][1];
    retval.data[2][0] = data[0][2];
    retval.data[2][1] = data[1][2];


    return retval;
  }

  // TESTED: WORKS
  matrix33 rot_x(REAL alpha) {
    
    matrix33 R;
    R.data[0][0] = 1;
    R.data[0][1] = R.data[0][2] = R.data[1][0] = R.data[2][0] = 0;
    R.data[1][1] = R.data[2][2] = cos(alpha);
    R.data[1][2] = -sin(alpha);
    R.data[2][1] = sin(alpha);

    return R;

  }


  // TESTED: WORKS
  matrix33 rot_y(REAL alpha) {

    matrix33 R;
    R.data[1][1] = 1;
    R.data[0][1] = R.data[1][0] = R.data[1][2] = R.data[2][1] = 0;
    R.data[0][0] = R.data[2][2] = cos(alpha);
    R.data[2][0] = -sin(alpha);
    R.data[0][2] = sin(alpha);

    return R;

  }

  // TESTED: WORKS
  matrix33 rot_z(REAL alpha) {

    matrix33 R;
    R.data[2][2] = 1;
    R.data[0][2] = R.data[1][2] = R.data[2][0] = R.data[2][1] = 0;
    R.data[0][0] = R.data[1][1] = cos(alpha);
    R.data[0][1] = -sin(alpha);
    R.data[1][0] = sin(alpha);

    return R;

  }

  matrix33& matrix33::operator=(const matrix33& a) {
    unsigned int i,j;
    for ( i = 0 ; i < 3 ; i ++ )
      for ( j = 0 ; j < 3 ; j++ )
	data[i][j] = a.data[i][j];
    return *this;
  }


  matrix33 operator+(const matrix33&a,const matrix33& b) {
    matrix33 r;
    unsigned int i,j;
    
    for ( i = 0 ; i < 3 ; i++ ) 
      for ( j = 0 ; j < 3 ; j++ ) 
	r.data[i][j]  = a.data[i][j] + b.data[i][j];


    return r;
  }

  s2_point operator*(const matrix33& m,const s2_point& v) {

    s2_point r;
    r.x = m.data[0][0]*v.x+ m.data[0][1]*v.y + m.data[0][2]*v.z;
    r.y = m.data[1][0]*v.x+ m.data[1][1]*v.y + m.data[1][2]*v.z;
    r.z = m.data[2][0]*v.x+ m.data[2][1]*v.y + m.data[2][2]*v.z;
    return r;

  }


  matrix::matrix(unsigned int R,unsigned int C) {
    r = R;
    c = C;
    allocate_memory();


  }

  void matrix::allocate_memory() { 
    unsigned int i;
    data = new REAL*[r];
    for ( i = 0 ; i < r ; i++ ) {
      data[i] = new REAL[c];
      fill(data[i],data[i]+c,0);
    }
  }

  void matrix::dispose_memory() {

    for ( unsigned int i = 0 ; i < r ; i++ )
      delete []data[i];
    delete []data;
    
  }
  
  matrix::matrix(const matrix33&a) {
    r = c = 3;
    allocate_memory();
    unsigned int i,j;
    for ( i = 0 ; i < r ; i++ )
      for ( j = 0 ; j < c ; j++ )
	data[i][j] = a.data[i][j];
  }

  void matrix::resize(unsigned int R,unsigned int C) {

    if ( ( r == R ) && ( c == C ) )
      return;
    dispose_memory();
    r = R;
    c = C;
    allocate_memory();
  }
  
  matrix::~matrix() {
    dispose_memory();
  }

  // TESTED: WORKS
  REAL matrix::norm() const {
    REAL res = 0; 
    unsigned int i,j;
    for (  i = 0 ; i < r ; i++ )
      for (  j = 0 ; j < c ; j++ )
	res += data[i][j]*data[i][j];
    res = sqrt(res);
    return res;
  }

  // TESTED: WORKS
  void matrix::normalize() {
    REAL res = norm(); 
    unsigned int i,j;
    for (  i = 0 ; i < r ; i++ )
      for (  j = 0 ; j < c ; j++ )
	data[i][j] /= res;
    
  }

  matrix& matrix::operator=(const matrix& s) {

    assert( (s.r == r ) && (s.c == c));
    for ( unsigned int i = 0 ; i < r ; i++ )
      for ( unsigned int j = 0 ; j < c ; j++ )
	data[i][j] = s.data[i][j];

    return *this;

  }

  // TESTED: WORKS
  matrix& matrix::operator*=(REAL f) {
    for ( unsigned int i = 0 ; i < r ; i++ )
      for ( unsigned int j = 0 ; j < c ; j++ )
	data[i][j] *= f;
    return *this;
  }


  matrix operator*(const matrix& A,const matrix&B) {

    assert(A.c == B.r);
    matrix r(A.r,B.c);
    unsigned int i,j,k;
    for ( i = 0 ; i < r.r ; i++ ) 
      for ( j = 0 ; j < r.c ; j++ )
	for ( k = 0 ; k < A.c ; k++ )
	  r.data[i][j] += A.data[i][k]*B.data[k][j];
    return r;
  }

  matrix operator-(const matrix&A,const matrix&B) {

    assert( (A.r == B.r ) && ( A.c == B.c ) );
    matrix r(A.r,A.c);
    unsigned int i,j;
    for ( i = 0 ; i < r.r ; i++ )
      for ( j =  0 ; j < r.c ; j++ ) 
	r.data[i][j] = A.data[i][j]-B.data[i][j];
    return r;
    
  }

  // TESTED: WORKS
  void matrix::save_to_file(const char*fname)const {
    unsigned int i,j;
    ofstream of(fname);
    for ( i = 0 ; i < r ; i++ ) {
      for ( j = 0 ; j < c ; j++ )
	of << data[i][j] << " " ;
      of << endl;
    }
    of.close();

  }

  // TESTED: WORKS
  matrix pointwise_prod(const matrix&a,const matrix&b) {

    assert((a.r == b.r ) && ( a.c == b.c) );
    matrix r(a.r,a.c);
    for ( unsigned int i = 0 ; i < a.r ; i++ )
      for ( unsigned int j = 0 ; j < a.c ; j++ )
	r.data[i][j] = a.data[i][j]*b.data[i][j];

    return r;

  }

  // TESTED: WORKS   
  matrix33 rot_axis_angle(const s2_point& axis,REAL angle) {
    REAL n = norm(axis);
    s2_point a = axis;
    a.x /= n; a.y /= n ; a.z /= n;

    matrix33 R,H;
    R.set_identity();
    H.data[0][0] = H.data[1][1] = H.data[2][2] = 0;
    H.data[0][1] = -a.z;
    H.data[0][2] = a.y;
    H.data[1][0] = a.z;
    H.data[1][2] = -a.x;
    H.data[2][0] = -a.y;
    H.data[2][1] = a.x;

    R = R + sinf(angle)*H + (1-cosf(angle))*H*H;

    return R;
    
  }

  s2_point normalize(const point& s) {
    
    s2_point v;
    REAL n = norm(s);
    v.x = s.x/n; v.y = s.y/n ; v.z = s.z/n;
    return v;

  }

  s2_point get_one_perpendicular(const s2_point& v) {
    s2_point r;
    s2_point a;

    if ( fabs(dot(v,s2_point(1,0,0))) < 0.9 ) {
      a = cross(v,s2_point(1,0,0));
      r = normalize(a);
    }
    else
      r = normalize(cross(v,s2_point(0,1,1)));
    
    assert(fabs(dot(r,v))< 0.00001);
    return r;

  }


  // TESTED: WORKS
  REAL dot_vector(REAL *a,REAL *b, int n) {

     int i;
    REAL r = 0;
    for ( i = 0 ; i < n ; i++ )
      r += a[i]*b[i];

    return r;

  }

  void rot_matrix_to_axis_angle(const matrix33&R,s2_point& v,REAL& theta) {
    
    theta = acosf( 0.5 * (R.data[0][0]+R.data[1][1]+R.data[2][2] - 1 ) );
    REAL m = 1/(2*sinf(theta));
    v.x = m * ( R.data[2][1] - R.data[1][2] );
    v.y = m * ( R.data[0][2] - R.data[2][0] );
    v.z = m * ( R.data[1][0] - R.data[0][1] );
    
  }

  // important: space for c should have been already allocated
  // TESTED: WORKS
  void circular_cross_correlation(REAL*a,REAL *b,REAL *c,int nel, int *lags, int nlags) {

    REAL *aa = new REAL[2*nel];
    REAL *aa_part = new REAL[nel];
    int k,j,i,lag;
    
    vector<REAL> v;

    for ( k = 0 ; k < nel ; k++ )
      aa[k] = aa[k+nel] = a[k];
    for ( k = 0 ; k < nlags ; k++ ) {
      lag = mod(lags[k],nel);
      for ( j = lag,i=0 ; j <= lag + nel - 1 ; j++ , i++ )
	aa_part[i] = aa[j];
     
      c[k] = dot_vector(aa_part,b,nel);
    }

    delete []aa_part;
    delete []aa;
  }

  // TESTED: WORKS
  REAL med3(REAL a,REAL b,REAL c) {

    if ( ( a < b ) && ( a < c ) )  {
      if ( b < c ) 
	return b;
      else return c;
    }
    if ( (  b < a ) && ( b < c ) )  {
      if ( a < c ) 
	return a;
      else
	return c;
    }
    if ( a < b ) 
      return a;
    else return b;

  }

  // TESTED: WORKS
  void medfilt1(const std::vector<REAL>&a,std::vector<REAL>& out,int order) {
    
    assert(order==3);

    int i,asize = a.size();
    assert(asize >= 3);
    out.resize(asize);
    
    out[0] = med3(0,a[0],a[1]);
    out[asize-1] = med3(0,a[asize-1],a[asize-2]);
    for ( i = 1 ; i < asize - 1 ; i++ )
      out[i] = med3(a[i-1],a[i],a[i+1]);

  }

  
  // Cross correlation between a and b is put in c. No normalization occurs.
  // TESTED: WORKS
  void cross_correlation(const vector<REAL>&a,const vector<REAL>&b,vector<REAL>& c,vector<int>& lags) {

    int s1 = a.size(), s2 = b.size();  
    int s = max(s1,s2);
    unsigned int outsize = s1 + s2 - 1;
    int tau,k;
    
    lags.resize(2*s-1);
    for ( k = 0 , tau = -s+1; k < 2*s -1 ; k++ ,tau++ )
      lags[k] = tau;
    
    c.resize(outsize);
    fill(c.begin(),c.end(),0);

    for ( tau = -s2 + 1 ; tau < s1 ; tau++ ) {
      for ( k = 0 ; k < s2 ; k++ ) 
	if ( ( k+tau >= 0 ) && ( k + tau < s1 )  )
	  c[tau+s2-1] += a[k+tau] * b[k];
    }
  

  }
  
  // TESTED: works
  void interp2_size(unsigned int r,unsigned int c,unsigned int order,
		    unsigned int& newrow,unsigned int& newcol) {
    newrow = r;
    newcol = c;
    for ( unsigned int i = 0 ; i < order ; i++ ) {
      newrow = 2 * newrow - 1;
      newcol = 2 * newcol - 1;
    }
  }

  // finds local maxima in a vector
  // if consider_extreme is true elements in 0 and n-1 may be maxima as well
  deque<corr_peak> find_maxima(const vector<REAL>& a,bool consider_extreme) {
    
    deque<corr_peak> retval;
    unsigned int size = a.size();
    corr_peak tmp;
    for ( unsigned int i = 1 ; i < size-1 ; i++ )
      if ( ( a[i]>=a[i-1] ) && ( a[i] >= a[i+1] ) ) {
	tmp.value = a[i];
	tmp.pos = i;
	retval.push_back(tmp);
      }

    if ( consider_extreme ) {
      if ( a[0] >= a[1] ) {
	tmp.value = a[0];
	tmp.pos = 0;
	retval.push_front(tmp);
      }
      if ( a[size-1] >= a[size-2] ) {
	tmp.value = a[size-1];
	tmp.pos = size-1;
	retval.push_back(tmp);
      }
    }

    return retval;

  }

  void solve_least_squares(const list<point>& ML,const list<REAL>&MY,const list<REAL>&weights,matrix& t_est,matrix33&cov,REAL&chi) {

    unsigned int i,n;
    n = weights.size();
    matrix W(n,n);
    matrix M(n,3);
    matrix MT(3,n);
    matrix Y(n,1);
    matrix residuals(n,1);
    list<REAL>::const_iterator lit;
    matrix cov_prod(3,3);
    
    for ( i = 0 , lit = weights.begin() ; i<n ; i++, lit++ ) 
      W.data[i][i] = *lit;

    for ( i = 0 , lit = MY.begin() ; i<n ; i++, lit++ ) 
      Y.data[i][0] = *lit;

    list<point>::const_iterator pit;
    for ( i = 0 , pit = ML.begin() ; i < n ; i++ , pit++ ) {
      M.data[i][0] = MT.data[0][i] = pit->x;
      M.data[i][1] = MT.data[1][i] = pit->y;
      M.data[i][2] = MT.data[2][i] = pit->z;
    }

    cov_prod = MT * W * M;
    cov = inv33(cov_prod);

    t_est.resize(3,1);
    t_est = matrix(cov) * MT * W * Y;
    residuals = M*t_est - Y;

    chi = 0;
    for ( i = 0 ; i < n ; i++ )
      chi += residuals.data[i][0] * W.data[i][i] * residuals.data[i][0];

  }

  // This is my own version of interp2
  // Important: enough space should have been ALREADY allocated for dst
  // TESTED: WORKS
  void interp2(REAL **src,unsigned int r,unsigned int c,unsigned int order,REAL **dst)  {

    unsigned int newrow,newcol;
    interp2_size(r,c,order,newrow,newcol);

    unsigned int step = 1;
    step <<= order;
    REAL increase;

    unsigned int ni,nj,i,j,k;
    
    // first copy carry-over data
    for ( j = 0 , nj  = 0 ; j < c ; j++ , nj += step )
      for ( i = 0 , ni = 0 ; i < r ; i++ , ni += step )
	dst[nj][ni] = src[j][i];

    // interpolate rows
    for ( ni = 0 , i = 0 ; ni < newrow ;  ni+= step , i++ ) {
      for ( nj = 0 , j = 0 ; nj < newcol-1 ; nj += step , j++ ) {
	increase  = ( src[i][j+1] - src[i][j] ) / step;
	for ( k = 1 ; k < step ; k++ )
	  dst[ni][nj+k] = dst[ni][nj+k-1] + increase;
      }
    }

    // interpolate columns
    for ( nj = 0 , j = 0 ; nj < newcol ;  nj+= step , j++ ) {
      for ( ni = 0 , i = 0 ; ni < newrow-1 ; ni += step , i++ ) {
	increase  = ( src[i+1][j] - src[i][j] ) / step;
	for ( k = 1 ; k < step ; k++ )
	  dst[ni+k][nj] = dst[ni+k-1][nj] + increase;
      }
    }

    REAL ri,ci;
    REAL REALstep = step;
    
    // interpolate everything else
    for ( ni = 0 ; ni < newrow -1 ; ni += step )
      for ( nj = 0 ; nj < newcol - 1 ; nj += step ) {
	for ( i = 1 ; i < step ; i++ )
	  for ( j =  1 ; j < step ; j++ ) {
	    ri =  i / REALstep;
	    ci =  j / REALstep;
	    dst[ni+i][nj+j] = dst[ni][nj] + (dst[ni][nj+step]-dst[ni][nj])*ci + 
	      (dst[ni+step][nj]-dst[ni][nj])*ri + (dst[ni+step][nj+step]+dst[ni][nj]-dst[ni][nj+step]-dst[ni+step][nj])*ri*ci;
	  }
      }
    
       

  }



  h3d_params::h3d_params() {

    max_norm = 8;
    rho_cell_size = 0.14899;
    angular_cell_size_deg = 3;
    upsample = 1;
    guess_rot_max_peaks1 = 2;
    guess_rot_max_peaks2 = 6  ;
    guess_rot_first_corr_num_hyp = 3;
    guess_rot_max_correction_deg = 6;
    guess_rot_num_compensations = 4;
    guess_rot_merge_hyp_threshold_deg = 5;
    closest_peaks_threshold = 15;
    cyl_amplitude_deg = 150;
    cyl_num_alpha_steps  = 15;
    cyl_num_beta_steps = 256;
    guess_tran_max_directions = 8;
    guess_tran_max_num_hyp = 4;
    guess_tran_two_birds_one_stone = 1;
    hausdorff_num_points = 200;
    hausdorff_plausible_threshold = 0.3;
    debug_compensation_show = true;
    debug_true_R = false;
    debug_true_T = false;
  }



}
