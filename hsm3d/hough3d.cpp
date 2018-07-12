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

#include <cmath>
#include <algorithm>
#include <map>
#include <cassert>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <deque>
#include <ctime>
#include <sys/time.h>


using namespace std;


namespace h3d {
 

  // TESTED: works
  bool operator<(const s2_peak& a,const s2_peak& b) {
    return a.value > b.value;  // this is correct ; it enforces reverse sorting
  }
  
  // TESTED: works
  bool operator<(const rotation_hypothesis&a,const rotation_hypothesis&b) {
    return a.value > b.value; // this is correct ; it enforces reverse sorting
  }
  
  // TESTED: works
  bool operator<(const translation_hypothesis&a,const translation_hypothesis&b) {
    return a.value > b.value; // this is correct ; it enforces reverse sorting
  }

  bool operator<(const ranked_solution&a,const ranked_solution&b) {
    return a.value > b.value;
  }

  bool operator<(const new_ranked_solution&a,const new_ranked_solution&b) {
    return a.distance < b.distance;
  }
  
  // TESTED: works
  bool operator<(const corr_peak& a,const corr_peak&b) {
    return a.value > b.value;
  }

  
  // TESTED: works
  h3d_struct::h3d_struct(REAL a,REAL r,REAL rmin,REAL rmax) {
    angular_cell_size_deg = a;
    rho_min = rmin;
    rho_max = rmax;

    rho_cell_size = r;

    rho_ncells = static_cast<unsigned int>(round((rho_max-rho_min)/r));
    // cube_ncells = static_cast<unsigned int>(ceil((M_PI/2)/deg2rad(angular_cell_size_deg)));
    cube_ncells = static_cast<unsigned int>(ceil( 90.0 / angular_cell_size_deg));

    unsigned int i,j,k;

    cube_point cp;
    buffer_point bp;
    s2_point normal;

    for ( bp.face = 1 ; bp.face <=6 ; bp.face++ )
      for ( bp.i = 1 ; bp.i <= cube_ncells ; bp.i++ )
        for ( bp.j = 1 ; bp.j <= cube_ncells ; bp.j++ ) {
          cp = coords_cell_to_avg_cube(bp,cube_ncells);
          normal = coords_cube_to_s2(cp);
          cube_normals.insert(make_pair(bp,normal));
        }
    

    // allocate structure cube_hs
    for ( i = 0 ; i < 6 ; i++ ) {
      cube_hs[i] = new REAL* [cube_ncells];
      for ( j = 0 ; j < cube_ncells ; j++ )
	      cube_hs[i][j] = new REAL[cube_ncells];
      
      }
    // fill it with 0s
    for ( i = 0 ; i < 6 ; i++ )
      for ( j = 0 ; j < cube_ncells ; j++ )
        fill(cube_hs[i][j],cube_hs[i][j]+cube_ncells,0);


    // allocate structure cube_ht
    for ( i = 0 ; i < 6 ; i++ ) {
      cube_ht[i] = new REAL** [cube_ncells];
      for ( j = 0 ; j < cube_ncells ; j++ ) {
	      cube_ht[i][j] = new REAL*[cube_ncells];
	      for ( k = 0 ; k < cube_ncells ; k++ )
	        cube_ht[i][j][k] = new REAL[rho_ncells];
      }
    
    }
    // fill it with 0s
    for ( i = 0 ; i < 6 ; i++ )
      for ( j = 0 ; j < cube_ncells ; j++ )
	     for ( k = 0 ; k < cube_ncells ; k++ )
	     fill(cube_ht[i][j][k],cube_ht[i][j][k]+rho_ncells,0);
  }

  h3d_struct::~h3d_struct() {

    unsigned int i,j,k;
    for ( i = 0 ; i < 6 ; i++ ) {
      for ( j = 0 ; j < cube_ncells ; j++ )
	      delete []cube_hs[i][j];
        delete []cube_hs[i];
    }

    for ( i = 0 ; i < 6 ; i++ ) {
      for ( j = 0 ; j < cube_ncells ; j++ ) {
        for ( k = 0 ; k < cube_ncells ; k++ ) 
	        delete []cube_ht[i][j][k];
	        delete []cube_ht[i][j];
      }
      delete []cube_ht[i];
    }
  }

  // TESTED: works
  void h3d_struct::compute_cube_ht(const vector<point>& c) {

    ptsImg = c;
    buffer_point bp;
    s2_point normal;

    map<buffer_point,s2_point>::iterator it;
    vector<point>::const_iterator pit;
    point p;
    REAL rhos,alpha;
    int r;
    int counter = 0;

    for ( bp.face = 1 ; bp.face <= 6 ; bp.face++ ) {
      for ( bp.i = 1 ; bp.i <= cube_ncells ; bp.i++ ) {
	for ( bp.j = 1 ; bp.j <= cube_ncells ; bp.j++ ) {
	  it = cube_normals.find(bp);
	  assert(it != cube_normals.end());
	  normal = it->second;
	  for ( pit = c.begin() ; pit != c.end() ; pit++ ) {
	    p = *pit;
	    rhos = p.x*normal.x + p.y * normal.y + p.z * normal.z;

	    rho2index(*this,rhos,r,alpha);
	
	      // cout << rhos << " " << r << " " << alpha << endl;
    
	    //  if (params.ht_points_are_from_scan) {
	      if  ( rhos <= EPSILON ) {
		counter++;
		continue;
	      }
		
	      //  }
	  
	    if ( r > 0 ) {
	      cube_ht[bp.face-1][bp.i-1][bp.j-1][r] += 1 - fabs(alpha);
	      if ( ( alpha > 0 ) && ( r < static_cast<int>(rho_ncells) - 1 ) )
		cube_ht[bp.face-1][bp.i-1][bp.j-1][r+1] += abs(alpha);
	      if ( ( alpha < 0 ) && ( r > 0 )  )
		cube_ht[bp.face-1][bp.i-1][bp.j-1][r-1] += abs(alpha);
	    }
	  }
	 
	}
      }
    }
   
    REAL sum;
    unsigned int k;

    for ( bp.face = 0 ; bp.face < 6 ; bp.face++ ) {
      for ( bp.i = 0 ; bp.i < cube_ncells ; bp.i++ )
	for ( bp.j = 0 ; bp.j < cube_ncells ; bp.j++ ) {
	  for ( sum = 0 , k=0 ; k < rho_ncells ; k++ )
	    sum+= cube_ht[bp.face][bp.i][bp.j][k]*cube_ht[bp.face][bp.i][bp.j][k];
	  cube_hs[bp.face][bp.i][bp.j] = sum;
	}
    }


  }

  void h3d_struct::compute_cube_ht_normals(const scan3D& scan,const matrix33& R, const point& T) {

    unsigned int ray,k;
    point p;
    s2_point alpha,direction;
    REAL rho,weight,a,sum;
    int r;
    buffer_point bp;
    map<buffer_point,s2_point>::iterator cube_it;

    for ( k = 0 ; k < scan.nscans ; k++ ) 
      for (  ray = 0 ; ray < scan.nrays ; ray++ ) {
	if ( ! scan.alpha_valid[k][ray] )
	  continue;

	p = scan.p3d[k][ray];
	alpha = scan.alpha[k][ray];

	p = R*p + T;
	alpha = R * alpha;

	bp = coords_s2_to_cell(alpha,cube_ncells);
	cube_it = cube_normals.find(bp);
	assert(cube_it != cube_normals.end());
	direction = cube_it->second;
	rho = dot(direction,p);
	rho2index(*this,rho,r,a);
	weight = scan.area[k][ray];

	if ( r > 0 ) {
	  cube_ht[bp.face-1][bp.i-1][bp.j-1][r] += (1 - fabs(a))*weight;
	  if ( ( a > 0 ) && ( static_cast<unsigned int>(r) < rho_ncells-1 ) )
	    cube_ht[bp.face-1][bp.i-1][bp.j-1][r+1] +=  fabs(a)*weight;
	  if ( ( a < 0 ) && ( r > 0 ) )
	    cube_ht[bp.face-1][bp.i-1][bp.j-1][r-1] +=  fabs(a)*weight;
	}
      }

    for ( bp.face = 0 ; bp.face < 6 ; bp.face++ ) {
      for ( bp.i = 0 ; bp.i < cube_ncells ; bp.i++ )
	for ( bp.j = 0 ; bp.j < cube_ncells ; bp.j++ ) {
	  for ( sum = 0 , k=0 ; k < rho_ncells ; k++ )
	    sum+= cube_ht[bp.face][bp.i][bp.j][k];//*cube_ht[bp.face][bp.i][bp.j][k];
	  cube_hs[bp.face][bp.i][bp.j] = sum;
	}
    }


  }


  // TESTED: works
  void h3d_struct::save_hs_to_file(const string& basename) {

    unsigned int face;
    ostringstream os;
    ofstream ofs;
    string fname;
    unsigned int i,j;
    for ( face = 0 ; face < 6 ; face++ ) {
      os << basename << face+1;
      fname = os.str();
      ofs.open(fname.c_str());
      os.str("");
      for ( i = 0 ; i < cube_ncells ; i++ ) {
	for ( j = 0 ; j < cube_ncells ; j++ )
	  ofs << cube_hs[face][i][j] << " ";
	ofs << endl;
      }
      ofs.close();
      
    }
  }
    
  
  // TESTED: WORKS
  h3d_struct h3d_struct::upsample(unsigned int ratio) const {

    unsigned int newrow,newcol,face;
    interp2_size(cube_ncells,cube_ncells,ratio,newrow,newcol);

    REAL new_angular_cell_size_deg = angular_cell_size_deg * cube_ncells / newrow;

    h3d_struct retval(new_angular_cell_size_deg,rho_cell_size,rho_min,rho_max);

    assert(retval.cube_ncells == newrow);

    for ( face = 0 ; face < 6 ; face++ ) 
      interp2(cube_hs[face],cube_ncells,cube_ncells,ratio,retval.cube_hs[face]);
     
    return retval;

  }


  // TESTED: WORKS
  list<buffer_peak> h3d_struct::matrix_local_maxima(REAL **m) {

    list<buffer_peak> retval;
    buffer_peak to_insert;
    unsigned int i,j;
    int val_count,k;
    REAL vals[8];
    bool is_peak;

    for ( i = 0 ; i < cube_ncells ; i++ )
      for ( j = 0 ; j < cube_ncells ; j++ ) {
	if ( m[i][j] == 0 ) 
	  continue;
	val_count = 0;
	if ( i > 0 ) vals[val_count++] = m[i-1][j];
	if ( j > 0 ) vals[val_count++] = m[i][j-1];
	if ( i < cube_ncells-1 ) vals[val_count++] = m[i+1][j];
	if ( j < cube_ncells-1 ) vals[val_count++] = m[i][j+1];
	
	if ( ( i > 0 ) && ( j > 0 ) )  vals[val_count++] = m[i-1][j-1];
	if ( ( i > 0 ) && ( j < cube_ncells-1 ) ) vals[val_count++] = m[i-1][j+1];
	if ( ( i < cube_ncells - 1 ) && ( j < cube_ncells - 1 ) ) vals[val_count++] = m[i+1][j+1];
	if ( ( i < cube_ncells - 1 ) && ( j > 0 ) ) vals[val_count++] = m[i+1][j-1];

	
	is_peak = true;
	for ( k = 0 ; k < val_count ; k++ ) {
	  if  ( vals[k] > m[i][j] ) { //|| ( REAL_compare(vals[k],m[i][j] ) ) ) {
	    is_peak = false;
	    break; // shortcut
	  }
	}
	if ( is_peak ) {
	  to_insert.p.face = 1; // useless -- just to get rid of a warning
	  to_insert.p.i = i+1;
	  to_insert.p.j = j+1;
	  to_insert.value = m[i][j];
	  retval.push_back(to_insert);
	}
      }

    return retval;

  }

  bool closer_than(const list<rotation_hypothesis>& Rs,const rotation_hypothesis& R,int t) {
    list<rotation_hypothesis>::const_iterator it;
    matrix33 T;
    s2_point v;
    REAL theta;
    for ( it = Rs.begin() ; it != Rs.end() ; it++ ) {
      T = R.rotation.transpose();
      T = it->rotation * T;
      rot_matrix_to_axis_angle(T,v,theta);
      if ( rad2deg(fabs(theta)) < t ) 
	return true;
    }
    return false;

  }
  
  bool closer_than_deg_or_parallel(const list<s2_point>& peaks,const s2_point& test,int max_deg) {

    list<s2_point>::const_iterator it;
    for ( it = peaks.begin() ; it != peaks.end() ; it++ ) {
      if ( angle_between_vectors(test,*it) < deg2rad(max_deg) )
	return true;
      if ( angle_between_vectors(-1*test,*it) < deg2rad(max_deg) )
	return true;
    }
    return false;

  }
  // TESTED: WORKS
  list<unsigned int> h3d_struct::filter_results(const vector<s2_peak>& v,REAL t) {
    list<unsigned int> retval;
    retval.push_back(0);
    unsigned int i;
    s2_peak vold,vnew; 
    vold = v[0];
    for ( i = 1 ; i < v.size() ; i++ )  {
      vnew = v[i];
      if ( apart(vold.p,vnew.p,t) ) {
	retval.push_back(i);
	vold = v[i];
      }
    }
    return retval;
  }

  // TESTED: WORKS
  list<s2_peak> h3d_struct::find_hs_peak() {
    list<buffer_peak> global_list,tmp_list;
    list<buffer_peak>::iterator pit;
    unsigned int face;
    for ( face = 0 ; face < 6 ; face++ ) {
      tmp_list.clear();
      tmp_list = matrix_local_maxima(cube_hs[face]);
      for ( pit = tmp_list.begin() ; pit != tmp_list.end() ; pit++ ) {
      	pit->p.face = face+1;
      	pit->value /= 10000; // actually useless; needed only for debug
      	global_list.push_back(*pit);
      }
    }

    vector<s2_peak> retval(global_list.size());
    s2_peak to_insert;
    cube_point tmp;
    unsigned int k = 0;
    for ( pit = global_list.begin() ; pit != global_list.end() ; pit++ ) {
      to_insert.value = pit->value;
      tmp = coords_cell_to_avg_cube(pit->p,cube_ncells);
      to_insert.p = coords_cube_to_s2(tmp);
      retval[k++] = to_insert;
    }

    stable_sort(retval.begin(),retval.end());  

    list<unsigned int> filtered = filter_results(retval,10);
    list<s2_peak> r;
    for (list<unsigned int>::iterator lit = filtered.begin() ; lit != filtered.end() ; lit++ )
      r.push_back(retval[*lit]);
    
    return r;
  }


  // TESTED: 
  // first element in  the couple is  theta
  list<pair<REAL,REAL> > h3d_struct::compensate_one_rotation(const h3d_struct&a,const s2_point&axis,const matrix33& current_guess,REAL max_correction,unsigned int spice) {

    
    s2_point support = get_one_perpendicular(axis);

    unsigned int num_alpha_steps = params.cyl_num_alpha_steps;
    unsigned int num_beta_steps = params.cyl_num_beta_steps;
    REAL amplitude = deg2rad(params.cyl_amplitude_deg);
    bool weight_by_alpha = true;
 
    matrix cyl1(num_alpha_steps,num_beta_steps);
    matrix cyl2(num_alpha_steps,num_beta_steps);

    cyl1 = project_hs_on_cylinder(*this,axis,support,num_alpha_steps,num_beta_steps,amplitude,weight_by_alpha,spice);
    cyl2 = project_hs_on_cylinder(a,current_guess*axis,current_guess*support,num_alpha_steps,num_beta_steps,amplitude,weight_by_alpha,spice);  


    cyl1.normalize();
    cyl2.normalize();


    unsigned int max_correction_step = ceil(max_correction * num_beta_steps/(2*M_PI));
    unsigned int nlags = 2*max_correction_step + 1;
    int *lags = new int[nlags];
    unsigned int i;
    int val;
    for ( i = 0 , val = -max_correction_step ; i < nlags ; i++ , val++)
      lags[i]=val;
    
    unsigned int alpha;
    REAL *corr = new REAL[nlags];
    REAL *tmp_corr = new REAL [nlags];
    fill(corr,corr+nlags,0);
    for ( alpha = 0 ; alpha < num_alpha_steps ; alpha++ ) {
      circular_cross_correlation(cyl2.data[alpha],cyl1.data[alpha],tmp_corr,num_beta_steps,lags,nlags);
      for ( i = 0 ; i < nlags ; i++ ) 
	corr[i] += tmp_corr[i];
      
    }

    vector<REAL> v(corr,corr+nlags);
    //save_vector(v,"save_correlation");

    deque<corr_peak> pos_max;
    /* corr_peak tmp;
    for ( unsigned int i = 1 ; i < nlags-1 ; i++ )
      if ( ( corr[i]>corr[i-1] ) && ( corr[i] > corr[i+1] ) ) {
	tmp.value = corr[i];
	tmp.pos = i;
	pos_max.push_back(tmp);
      }
    */
    pos_max = find_maxima(v,false);
    sort(pos_max.begin(),pos_max.end());

    list<pair<REAL,REAL> >  retval;
    for ( deque<corr_peak>::iterator it = pos_max.begin() ; it != pos_max.end() ; it++ ) 
      retval.push_back(make_pair(lags[it->pos]*2*M_PI/num_beta_steps,it->value));
    


    delete []corr;
    delete []tmp_corr;

    delete []lags;
    return retval;

  }

  

  // TESTED: OK
  list<rotation_hypothesis> h3d_struct::guess_rotation(const h3d_struct& ht2) {

    deque<rotation_hypothesis> hypothesis;
  

    h3d_struct ht1u = upsample(params.upsample);
    h3d_struct ht2u = ht2.upsample(params.upsample);

    list<s2_peak> p1u = ht1u.find_hs_peak();
    list<s2_peak> p2u = ht2u.find_hs_peak();
    deque<s2_point> axes;

    unsigned int n1,n2,peak1,peak2,i,spice,a;
    
    list<s2_peak>::iterator p1,p2;

    s2_point m1,m2,cross_m1m2,tmp,axis;
    matrix33 r1,r2,Rest,r_i;
    REAL max_correction;
    rotation_hypothesis tmprot;

    n1 = min(static_cast<unsigned int>(p1u.size()),params.guess_rot_max_peaks1);
    n2 = min(static_cast<unsigned int>(p1u.size()),params.guess_rot_max_peaks2);

    list<pair<REAL,REAL> > comp1,comp2;
    list<pair<REAL,REAL> >::iterator c_it;

    for ( p1 = p1u.begin() , peak1 = 0; peak1 < n1 ; p1++ , peak1++)
      for ( p2 = p2u.begin() , peak2 = 0 ; peak2 < n2 ; p2++ , peak2++) {
	
      	m1 = p1->p;
      	m2 = p2->p;

	//special cases
	if ( dot(m1,m2) > 0.9999 )
	  r1.set_identity();
	else if ( dot(m1,m2) < -0.9999 )  // should cache it...
	  r1 = rot_axis_angle(get_one_perpendicular(m1),M_PI);
	else
	  r1 = rot_axis_angle(cross(m1,m2),angle_between_vectors(m1,m2));

	max_correction = static_cast<REAL>(M_PI);

	comp1 = ht1u.compensate_one_rotation(ht2u,m1,r1,max_correction,0);
	
	unsigned int howmany = min(static_cast<unsigned int>(comp1.size()),params.guess_rot_first_corr_num_hyp);
	for (  a = 0 , c_it = comp1.begin () ; a < howmany ; a++ , c_it++ ) {
	  r2 = rot_axis_angle(r1*m1,c_it->first);
	  Rest = r2 * r1;
	
	  axes.clear();
	  axes.push_back(get_one_perpendicular(m1));
	  axes.push_back(cross(m1,axes.front()));
	  axes.push_back(m1);
	  //  axes.push_back(axes.front());
		       

	  for ( i = axes.size() + 1; i <= params.guess_rot_num_compensations ; i++ ) {
	    tmp = rot_axis_angle(axes[i-3],deg2rad(30))*normalize(cross(axes[i-2],axes[i-3]));
	    //  tmp = normalize(cross(axes[i-2],axes[i-3]));
	    axes.push_back(tmp);
	  }

	  for ( i = 0 ; i < params.guess_rot_num_compensations ; i++ ) {
	    axis = axes[i];
	    max_correction = deg2rad(params.guess_rot_max_correction_deg);
	    spice = i+1;
	    comp2 = ht1u.compensate_one_rotation(ht2u,axis,Rest,max_correction,spice);

	    if ( comp2.size() > 0 )  {
	      r_i = rot_axis_angle(Rest*axis,comp2.front().first);
	      Rest = r_i * Rest;
	    }
	  }

	  tmprot.rotation = Rest;
	  tmprot.value = c_it->second;
	  hypothesis.push_back(tmprot);
	}
      }
  
    stable_sort(hypothesis.begin(),hypothesis.end());
	
      
    list<rotation_hypothesis> retval;
    
    unsigned int num_hyp = hypothesis.size();
    for ( i = 0 ; i < num_hyp ; i++ ) {
      if ( ! closer_than(retval,hypothesis[i],params.guess_rot_merge_hyp_threshold_deg) ) 
	retval.push_back(hypothesis[i]);
    }
    
    return retval;
    
  }

  // TESTED:
  deque<translation_hypothesis> h3d_struct::guess_translation(const h3d_struct& ht2,const matrix33& Rest) {

    list<s2_peak> peaks = find_hs_peak();
    list<s2_point> peaks_used;
    list<s2_peak>::iterator pit;
    s2_point peak;
    s2_point other_direction, other_direction_wished;
    buffer_point bp1,bp2,bp1b,bp2b;
    vector<REAL> s1(rho_ncells),s2(ht2.rho_ncells),s1b(rho_ncells),s2b(ht2.rho_ncells);
    vector<REAL> xc_unsmoothed,xc,maxima_delta;
    vector<int> lags;
    int j,pcount;
    unsigned int num_dir,num_hyp,i,h;
    REAL mismatch;
    deque<corr_peak> maxima;
    deque<corr_peak>::iterator max_it;
    map<buffer_point,s2_point>::const_iterator m_it;
    deque<deque<obs_obj> > obs;
    deque<obs_obj> to_insert;
    obs_obj tmp_obj;

    for ( pit = peaks.begin() , pcount = 0 ; pit != peaks.end() ; pit++ , pcount++ ){
      peak = pit->p;
      if ( closer_than_deg_or_parallel(peaks_used,peak,params.closest_peaks_threshold) )
	continue;
      bp1 = coords_s2_to_cell(peak,cube_ncells);
      bp2 = coords_s2_to_cell(Rest*peak,ht2.cube_ncells);
      if ( ht2.cube_hs[bp2.face-1][bp2.i-1][bp2.j-1] <= 0.0001 * cube_hs[bp1.face-1][bp1.i-1][bp1.j-1])
	continue;

      peaks_used.push_back(peak);

      copy(cube_ht[bp1.face-1][bp1.i-1][bp1.j-1],cube_ht[bp1.face-1][bp1.i-1][bp1.j-1]+rho_ncells,s1.begin());
      copy(ht2.cube_ht[bp2.face-1][bp2.i-1][bp2.j-1],ht2.cube_ht[bp2.face-1][bp2.i-1][bp2.j-1]+ht2.rho_ncells,s2.begin());

      if ( params.guess_tran_two_birds_one_stone ) {
	bp1b = coords_s2_to_cell((-1)*peak,cube_ncells);
	bp2b = coords_s2_to_cell(Rest * ((-1)*peak),ht2.cube_ncells);

	copy(cube_ht[bp1b.face-1][bp1b.i-1][bp1b.j-1],cube_ht[bp1b.face-1][bp1b.i-1][bp1b.j-1]+rho_ncells,s1b.begin());
	copy(ht2.cube_ht[bp2b.face-1][bp2b.i-1][bp2b.j-1],ht2.cube_ht[bp2b.face-1][bp2b.i-1][bp2b.j-1]+ht2.rho_ncells,s2b.begin());

	reverse(s1b.begin(),s1b.end());
	reverse(s2b.begin(),s2b.end());
	
	for ( i = 0 ; i < rho_ncells ; i++ )
	  s1[i] = s1[i] - s1b[i];

	for ( i = 0 ; i < ht2.rho_ncells ; i++ )
	  s2[i] = s2[i] - s2b[i];


	m_it = ht2.cube_normals.find(bp2);
	assert(m_it != ht2.cube_normals.end());
	other_direction = m_it->second;
	other_direction_wished = Rest * peak;
	mismatch = angle_between_vectors(other_direction,other_direction_wished);
	

      }

      cross_correlation(s2,s1,xc_unsmoothed,lags);
      //save_vector(xc_unsmoothed,"../matlab/saved_cross_correlation");
      medfilt1(xc_unsmoothed,xc,3);
      //save_vector(xc,"../matlab/saved_cross_smoothed");
      
      maxima = find_maxima(xc,true);
      if ( maxima.size() == 0 ) 
	continue;
      stable_sort(maxima.begin(),maxima.end());

      save_deque(maxima,"../matlab/saved_maxima");
	
      maxima_delta.resize(maxima.size());
      i = 0;
      for ( max_it = maxima.begin()  ; max_it != maxima.end() ; max_it++ , i++ )
	maxima_delta[i] = lags[max_it->pos] * rho_cell_size;
    
      //save_vector(maxima_delta,"../matlab/saved_maxima_delta");

      num_dir = obs.size()+1;
      to_insert.clear();
      
      bool too_close;
      int prev;

      for ( j = 0 ; j < (static_cast<int> (maxima_delta.size())) ; j++ ) {
	
	if ( maxima[j].value <= 0.01 * maxima[0].value ) 
	  continue;
	
	too_close = false;
	for ( prev = 0 ; prev <= j-1 ; prev++ ) {
	  if ( fabs(maxima_delta[prev]-maxima_delta[j]) < 0.2 ) 
	    too_close = true;
	}
	if ( too_close ) 
	  continue;
	
	num_hyp = to_insert.size()+1;
	if ( num_hyp > params.guess_tran_max_num_hyp ) 
	  break;
	
	tmp_obj.M_row = Rest*peak;
	tmp_obj.Z_row = maxima_delta[j];
	tmp_obj.weight = maxima[j].value;
	
	to_insert.push_back(tmp_obj);
	
      }
      if ( to_insert.size() > 0 )
	obs.push_back(to_insert);

      if ( peaks_used.size() >= params.guess_tran_max_directions ) 
	break;
      
    }

    deque<deque<int> > index_combinations;
    deque<deque<int> > new_combinations;	
    deque<int> tmp_deque;
    index_combinations.push_back(tmp_deque);
    unsigned int ndirections = obs.size();
    for ( i = 0 ; i < ndirections ; i++ ) {
      num_hyp = obs[i].size();
      new_combinations.clear();
      for ( h = 1 ; h <= num_hyp ; h++ ) 
	for ( j = 1 ; j <= static_cast<int> (index_combinations.size()) ; j++ ) {
	  new_combinations.push_back(index_combinations[j-1]);
	  new_combinations.back().push_back(h-1);
	}
      index_combinations = new_combinations;
    }
    /*
    for ( i = 0 ; i < index_combinations.size() ; i++ )  {
      for ( int pp = 0 ; pp < index_combinations[i].size() ; pp++ ) 
	cout << index_combinations[i][pp] << "  ";
      cout << endl;
      }
    */

    
    deque<translation_hypothesis> retval;
    translation_hypothesis tmp;
    unsigned int ncombinations = index_combinations.size();
    deque<int> combination;
    list<point> M;
    list<REAL> weights,Z;
    list<REAL>::iterator W_it;
    matrix t_est(3,1);
    matrix33 cov;
    REAL chi,sum_weights;
    for ( unsigned int a = 0 ; a < ncombinations ; a++ ) {
      combination = index_combinations[a];
      M.clear();
      weights.clear();
      Z.clear();
      for ( i = 0 ; i < ndirections ; i++ ) {
	h = combination[i];
	if ( obs[i][h].weight > 0 ) {
	  M.push_back(obs[i][h].M_row);
	  Z.push_back(obs[i][h].Z_row);
	  weights.push_back(obs[i][h].weight);
	}
      }
      sum_weights = 0;
      for ( W_it = weights.begin() ; W_it != weights.end() ; W_it++ ) 
	sum_weights += *W_it;
      solve_least_squares(M,Z,weights,t_est,cov,chi);
      tmp.value = sum_weights;
      tmp.direction.x = t_est.data[0][0];
      tmp.direction.y = t_est.data[1][0];
      tmp.direction.z = t_est.data[2][0];
      retval.push_back(tmp);
    }
  
    sort(retval.begin(),retval.end());
  
    return retval;
  }

  // TESTE: WORKS
  view_score h3d_struct::scan_whats_your_view(const scan3D& scan,const point& p) {

    tt_struct s = scan_dir2tt(p);
    view_score retval;

    retval.singular = s.singular;

    REAL *thetas = new REAL[scan.nrays];
    REAL *tilts = new REAL[scan.nscans];
    copy(scan.thetas_rad,scan.thetas_rad+scan.nrays,thetas);
    copy(scan.tilt_rad,scan.tilt_rad+scan.nscans,tilts);
    /*
    unsigned int i;
    for ( i = 0 ; i < scan.nrays ; i++ )
      thetas[i] = deg2rad(scan.thetas_deg[i]);
    for ( i = 0 ; i < scan.nscans ; i++ ) 
      tilts[i] = deg2rad(scan.tilt_deg[i]);
    */
    REAL epsilon = 0.0001;
    thetas[0] -= epsilon;
    thetas[scan.nrays-1] += epsilon;
    
    tilts[0] -= epsilon;
    tilts[scan.nscans-1] += epsilon;

    if ( ( s.theta < thetas[0] ) || 
	 ( s.theta > thetas[scan.nrays-1] ) || 
	 ( s.tilt < tilts[0] )  || 
	 ( s.tilt > tilts[scan.nscans-1] ) ) {
      retval.visible = false;
      retval.rho = -50000;
      delete []thetas;
      delete []tilts;
      return retval;
    } 

    unsigned int index_theta =  1 + round( (scan.nrays-1)*(s.theta-thetas[0])/(thetas[scan.nrays-1]-thetas[0]));
    unsigned int index_tilt = 1 + round( (scan.nscans-1)*(s.tilt-tilts[0])/(tilts[scan.nscans-1]-tilts[0]));

    assert( ( index_theta >= 1 ) && ( index_theta <= scan.nrays ) );
    assert( ( index_tilt >= 1 ) && ( index_tilt <= scan.nscans ) );
    
    retval.rho = scan.readings[index_tilt-1][index_theta-1];
    retval.visible = true;

    delete []thetas;
    delete []tilts;

    return retval;

  }
  

  vector<point> h3d_struct::sample_points(const scan3D&s,unsigned int npoints) {

    REAL tot_area=0;
    unsigned int i,k,j;
    unsigned int p,n = s.nscans*s.nrays,r,l,m,pos;
    REAL *probs = new REAL[n];
    bool *sampled = new bool[n];
    REAL r_num;

    fill(sampled,sampled+n,false);

    for ( k = 0 ; k < s.nscans ; k++ )
      for ( i = 0 ; i < s.nrays ; i++ )
	      tot_area += s.area[k][i];

    j = 0;
    for ( k = 0 ; k < s.nscans ; k++ )
      for ( i = 0 ; i < s.nrays ; i++ )
	      probs[j++] = s.area[k][i] / tot_area;

    for ( j = 1 ; j < n ; j++ )
      probs[j] += probs[j-1];
      
    probs[n-1] = 1;

    
    vector<point> v(npoints);
    
    timeval tv;
    gettimeofday(&tv,NULL);
    srandom(tv.tv_usec);

    unsigned int generated = 0;
    pos = 0;
    // roulette selection...
    for ( p = 0  ; p < npoints ; p++  ) {
      r_num = static_cast<REAL>(random())/RAND_MAX;
      if ( r_num <= probs[0] )
	      pos = 0;
      else {
	// binary search
      	r = n-1 ;
      	l = 0;
	      while ( l<=r ) {
	        m = (r+l)/2;
	        if ( ( probs[m-1] < r_num ) && ( probs[m] >= r_num ) ) {
	          pos = m;
	          l = r+1;
	        }
	        else if ( probs[m] < r_num ) 
	          l = m+1;
	        else 
	          r = m-1;
	      }
      }
      if ( ! sampled[pos] ) {  // avoid resampling the same point
	      v[generated] = s.points[pos];
	      sampled[pos] = true;
	      generated++;
      }
    }
    
    v.resize(generated);

    
    delete []probs;
    delete []sampled;
    return v;

  }
  
  vector<point> h3d_struct::choose_sample_points(const scan3D&s,unsigned int npoints){
    std::vector<point> v;
    return v;

  }
  
  vector<ranked_solution> h3d_struct::evaluate_solutions_hausdorff(const scan3D& scan1,const scan3D& scan2,const matrix33& add_R,const point& add_T,const list<roto_translation_hypothesis>& solutions) {

    vector<point> sampled_points;
    matrix33 R;
    point p,pw,T;
    unsigned int a,index;
    view_score s;
    REAL weight;
    vector<ranked_solution> retval(solutions.size());
    unsigned int npoints;
  
    sampled_points  = sample_points(scan2,params.hausdorff_num_points);
    npoints = sampled_points.size();
  
    list<roto_translation_hypothesis>::const_iterator s_it;
    for ( s_it = solutions.begin() , index = 0 ; s_it != solutions.end() ; s_it++ , index++) {

      R = s_it->rot.rotation;
      T = s_it->tran.direction;
      R = R.transpose();

      retval[index].rot = s_it->rot;
      retval[index].tran = s_it->tran;
      
      retval[index].num_out_of_field = 0;
      retval[index].num_behind_something = 0;
      retval[index].num_inconsistent = 0;
      retval[index].num_plausible = 0;
      for ( a = 0 ; a < npoints ; a++ ) {
      	p = sampled_points[a];
      	
      	pw = R * (p-T);
      	pw = add_R * pw + add_T;

      	s = scan_whats_your_view(scan1,pw);
      	weight = sqrt(dot(pw,pw));
      	if ( ! s.visible ) {
      	  retval[index].num_out_of_field++;
      	  continue;
      	}
      	if ( fabs(s.rho - weight) < params.hausdorff_plausible_threshold )
      	  retval[index].num_plausible += weight;
      	else if ( s.rho < weight ) 
      	  retval[index].num_behind_something += weight;
      	else
      	  retval[index].num_inconsistent += weight;
      }
      retval[index].value = retval[index].num_plausible;
    }
    
    sort(retval.begin(),retval.end());
    return retval;

  }

  vector<new_ranked_solution> h3d_struct::evaluate_solutions_hausdorff_pts(const vector<point>& c,const list<roto_translation_hypothesis>& solutions, ANNpointArray &kd_pts, int &npts){

    vector<new_ranked_solution> retval;
    new_ranked_solution it_t;
    point pit;
    matrix33 R;
    point T;
    int dim =3;
    ANNkd_tree* kdTree = new ANNkd_tree(kd_pts, npts, dim);
    
    list<roto_translation_hypothesis>::const_iterator it;
    // // Brute Force
    // for(it = solutions.begin();it!=solutions.end();++it){
    //   REAL min;
    //   REAL max = 0.0;
      
    //   R = it->rot.rotation;
    //   T = it->tran.direction;
    //   R=R.transpose();
    //   it_t.rot.rotation=R;
    //   it_t.tran.direction=T;
    //   bool flag = false;
    //   for(int i=0;i<c.size();++i){
          
    //     pit = R*(c[i]-T);
    //     for(int k=0;k<ptsImg.size();++k){
    //       REAL distance=sqrt(pow((ptsImg[k].x-pit.x),2)+pow((ptsImg[k].y-pit.y),2)+pow((ptsImg[k].z-pit.z),2));
    //       if(!flag) {
    //         min = distance; 
    //         flag=true;
    //       }
    //       if (distance<min) {
    //         min=distance;
    //       }  
           
    //     }
    //        if(min>max) max=min;
    //   }
    //   it_t.distance = min;
    //   retval.push_back(it_t);
    // }

    // KDTree Implementation
    for(it = solutions.begin();it!=solutions.end();++it){
      REAL max = 0;
      ANNpoint queryPt;
      ANNidxArray nnIdx;
      ANNdistArray dists;
      nnIdx = new ANNidx[1];
      dists = new ANNdist[1];
      R = it->rot.rotation;
      T = it->tran.direction;
      R = R.transpose();
      it_t.rot.rotation=R;
      it_t.tran.direction=T;
      for(int i=0;i<c.size();++i){
        pit = R*(c[i]-T);
        
        queryPt = annAllocPt(dim);
        queryPt[0] = pit.x;
        queryPt[1] = pit.y;
        queryPt[2] = pit.z;
        // cout << "(" << queryPt[0];
        // for (int l = 1; l < dim; l++) {
        //   cout << ", " << queryPt[l];
        // }
        // cout << ")\n";
        kdTree->annkSearch(queryPt,1,nnIdx,dists,0);
        dists[0] = sqrt(dists[0]);      // unsquare distance
        // cout << "\tNN:\tIndex\tDistance\n";
        // cout << "\t" << i << "\t" << nnIdx[0] << "\t" << dists[0] << "\n";
        // cout<<"Nearest Neighbor:\n";
        // cout << "(" << kd_pts[nnIdx[0]][0];
        // for (int m = 1; m < dim; m++) {
        //   cout << ", " << kd_pts[nnIdx[0]][m];
        // }
        // cout << ")\n";
        if(dists[0]>max) max = dists[0];

      }
      it_t.distance = max;
      retval.push_back(it_t);

      
      delete [] nnIdx;
      delete [] dists;
    }
    
    
    sort(retval.begin(),retval.end());
    return retval;

  }

}

