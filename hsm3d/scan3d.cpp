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


#include "scan3d.h"
#include "io.h"

#include <json/json.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

namespace h3d {

  // TESTED: WORKS
  scan3D::scan3D() {
    tilt_deg = thetas_deg = NULL;
    tilt_rad = thetas_rad = NULL;
    readings = area = NULL;
    points = NULL;
    nscans = nrays = 0;
    p3d = NULL;
    valid = NULL;
    alpha = NULL;
    alpha_valid = NULL;
  }

  // TESTED: WORKS
  scan3D::scan3D(unsigned int s,unsigned int r) {

    nscans = s;
    nrays = r;
    allocate_memory();

  }

  // This is to be used for benchmarking purposes only
  int scan3D::load_absolute_points_from_file(const char * fname) {
    
    if ( ( nscans == 0 ) && ( nrays == 0  ) )
      return 4;


    ifstream  ifs(fname);
    if ( ! ifs )
      return 1;

    unsigned int i,k;
    point p;
    
    //  tt_struct s;

    ifstream thetasfile("thetas_deg.txt");
    if ( ! thetasfile )
      return 1;
    for ( i = 0 ; i < nrays ; i++ )
      thetasfile >> thetas_deg[i];
    thetasfile.close();

    ifstream tiltsfile("tilts_deg.txt");
    if ( ! tiltsfile )
      return 1;
    for ( i = 0 ; i < nscans ; i++ )
      tiltsfile >> tilt_deg[i];
    tiltsfile.close();
   
    for ( i = 0 ; i < nscans ; i++ ) {
      for ( k = 0 ; k < nrays ; k++ ) { 
	  ifs >> p.x;
	  ifs >> p.y;
	  ifs >> p.z;
	  readings[i][k] = sqrt(dot(p,p));
	}
      }
  

    ifs.close();

    scan_convert_points();
    scan_find_normals();

    return 0;

  }
  int scan3D::load_points_from_file(const char* file){}


  // TESTED: WORKS
  void scan3D::allocate_memory() {
    unsigned int i;
    tilt_deg = new REAL[nscans];
    thetas_deg = new REAL[nrays];
    thetas_rad = new REAL[nrays];
    tilt_rad = new REAL[nscans];
    readings = new REAL*[nscans];
    area = new REAL*[nscans];
    p3d = new point*[nscans];
    valid = new bool*[nscans];
    alpha_valid = new bool*[nscans];
    alpha = new s2_point*[nscans];
    for ( i = 0 ; i < nscans ; i++ ) {
      readings[i] = new REAL[nrays];
      area[i] = new REAL[nrays];
      p3d[i] = new point[nrays];
      valid[i] = new bool[nrays];
      fill(valid[i],valid[i]+nrays,true);
      alpha_valid[i] = new bool[nrays];
      alpha[i] = new s2_point[nrays];
    }
    points = new point[nscans * nrays];
  }

  // TESTED:WORKS
  scan3D::~scan3D() {

    if ( tilt_deg )
      delete []tilt_deg;
    if ( thetas_deg )
      delete []thetas_deg;
    if ( thetas_rad )
      delete []thetas_rad;
    if ( tilt_rad )
      delete []tilt_rad;
    if ( nscans > 0 ) {
      for ( unsigned int i = 0 ; i < nscans ; i++ ) {
	if ( readings ) delete []readings[i];
	if ( area ) delete []area[i];
	if ( p3d ) delete []p3d[i];
	if ( valid ) delete []valid[i];
	if ( alpha_valid ) delete []alpha_valid[i];
	if ( alpha ) delete []alpha[i];
      }
    }
    if ( area ) delete []area;
    if ( readings ) delete []readings;
    if ( points  ) delete []points;
    if ( p3d ) delete []p3d;
    if ( valid ) delete []valid;
    if ( alpha_valid ) delete []alpha_valid;
    if ( alpha ) delete []alpha;
  }

  // return code:
  // 0: ok
  // -1: problem loading file
  // -2: file structure inconsistent with scan size
  // TESTED: WORKS
  int scan3D::load_scan_jsonfile(char *fname) {

    unsigned int i,j;
    json_object * jobj = json_object_from_file(fname);

    if ( ! jobj ) 
      return -1;

    json_object * thetas_deg_json = json_object_object_get(jobj,"thetas_deg");
    json_object * tilt_deg_json = json_object_object_get(jobj,"tilt_deg");
    json_object * readings_json = json_object_object_get(jobj,"readings");

    // error checking

    if ( ( ! thetas_deg_json ) || ( ! tilt_deg_json ) || ( ! readings_json ) )
      return -1;

    if ( json_object_get_type(thetas_deg_json) != json_type_array )
      return -1;

    if ( json_object_get_type(tilt_deg_json) != json_type_array )
      return -1;
    
    if ( json_object_get_type(readings_json) != json_type_array )
      return -1;
    

    if ( ( nrays == 0 ) && ( nscans == 0 ) ) {
      nrays = static_cast<unsigned int>(json_object_array_length(thetas_deg_json));
      nscans = static_cast<unsigned int>(json_object_array_length(tilt_deg_json));
      allocate_memory();
    }
    else {
      if ( ( nrays != static_cast<unsigned int>( json_object_array_length(thetas_deg_json) ) ) &&
	   ( nscans != static_cast<unsigned int>( json_object_array_length(tilt_deg_json) ) ) )
	return -2;
    }
    
    json_object *tmpjson;

    for ( i = 0 ; i < nscans ; i++ ) {
      tmpjson = json_object_array_get_idx(tilt_deg_json,i);
      tilt_deg[i] = (REAL)json_object_get_double(tmpjson);
      json_object_put(tmpjson);
    }

    for ( i = 0 ; i < nrays ; i++ ) {
      tmpjson = json_object_array_get_idx(thetas_deg_json,i);
      thetas_deg[i] = (REAL)json_object_get_double(tmpjson);
      json_object_put(tmpjson);
    }
    
    json_object * mrow;
  

    for ( i = 0 ; i < nscans ; i++ ) {
      mrow = json_object_array_get_idx(readings_json,i);
      assert( json_object_array_length(mrow)== static_cast<int>(nrays));
      for( j = 0 ; j < nrays ; j++ ) {
	tmpjson = json_object_array_get_idx(mrow,j);
	readings[i][j] = (REAL)json_object_get_double(tmpjson);
	json_object_put(tmpjson);
      }
      json_object_put(mrow);
  
    }
    
    json_object_put(thetas_deg_json);
    json_object_put(tilt_deg_json);
    json_object_put(readings_json);
    
    json_object_put(jobj);  
    
    scan_convert_points();
    scan_find_normals();
    
    return 0;

  }
  
  // TESTED: WORKS
  void scan3D::scan_convert_points(void) {

    if ( ( nrays == 0 )  || ( nscans == 0 ) )
      return;

    unsigned int i,k,np;
    REAL rho,theta,tilt;
    s2_point direction;
    point p_p3d;

    for ( k = 0 ; k < nscans ; k++ )
      tilt_rad[k] = deg2rad(tilt_deg[k]);
    for ( i = 0 ; i < nrays ; i++ )
      thetas_rad[i] = deg2rad(thetas_deg[i]);

    np = 0;
    for ( k = 0 ; k < nscans ; k++ )
      for ( i = 0 ; i < nrays ; i++ ) {
	rho = readings[k][i];
	if  (  rho > 0  ) {
	  theta = thetas_rad[i];
	  tilt = tilt_rad[k];

	  direction = (rot_y(tilt) * rot_z(theta)) * point(1,0,0);
	  p_p3d = rho * direction;

	  points[np++] = p_p3d;
	  p3d[k][i] = p_p3d;
	  valid[k][i] = true;
	}
	else {
	  p_p3d = point(0,0,0);
	  points[np++] = p_p3d;
	  p3d[k][i] = p_p3d;
	  valid[k][i] = false;
	}
	
      }

  }

  // TESTED: WORKS
  void scan3D::scan_find_normals() {

    unsigned int k,i;
    REAL rho,side1,side2;
    point v1,v2;
    s2_point n;

    for ( k = 0 ; k < nscans ; k++ ) {
      fill(area[k],area[k]+nrays,0);
      fill(alpha_valid[k],alpha_valid[k]+nrays,false);
    }

    for ( k = 1 ; k < nscans- 1 ; k++ )
      for ( i = 1 ; i < nrays - 1 ; i++ ) {
	
	rho = readings[k][i];
	if ( rho <= 0 )
	  continue;
	if ( ! ( valid[k+1][i] && valid[k-1][i] )  ) 
	  side1 = 0;
	else {
	  v1 = p3d[k+1][i] - p3d[k-1][i];
	  side1 = norm(v1);
	}
	if ( ! ( valid[k][i+1] && valid[k][i-1] )  ) 
	  side2 = 0;
	else {
	  v2 = p3d[k][i+1] - p3d[k][i-1];
	  side2 = norm(v2);
	}
	
	if ( ( side2 > 0 ) && ( side1 > 0 ) ) {
	  n = normalize(cross(v1,v2));
	  if ( dot(n,p3d[k][i]) > 0 )
	    n = -1.0 * n;
	  alpha[k][i] = n;
	  area[k][i] = side1*side2;
	  alpha_valid[k][i] = true;
	}
      }
    
  }
  

  // TESTED: WORKS
  void scan3D::save_scan_to_file(const char*fname) {

    string tiltname(fname),degname(fname),readingsname(fname);
    string areaname(fname),alphaname(fname),alphavalidname(fname);
    string pointsname(fname),p3dname(fname),validname(fname);

    tiltname = tiltname + "_tilt";
    degname = degname + "_deg";
    readingsname = readingsname + "_readings";
    areaname = areaname + "_area";
    alphaname = alphaname + "_alpha";
    alphavalidname = alphavalidname + "_alphavalid";
    pointsname = pointsname + "_points";
    p3dname = p3dname + "_p3d";
    validname = validname + "_valid";

    vector<REAL> tilt(tilt_deg,tilt_deg+nscans);
    save_vector(tilt,tiltname.c_str());
    vector<REAL> thetas(thetas_deg,thetas_deg+nrays);
    save_vector(thetas,degname.c_str());
    vector<point> pointsvec(points,points+(nscans*nrays));
    save_vector(pointsvec,pointsname.c_str());

    ofstream ofsreadings(readingsname.c_str());
    ofstream ofsarea(areaname.c_str());
    ofstream ofsalpha(alphaname.c_str());
    ofstream ofsalphavalid(alphavalidname.c_str());
    ofstream ofsp3d(p3dname.c_str());
    ofstream ofsvalid(validname.c_str());

    for (unsigned int i = 0 ; i < nscans ; i++ ) {
      for ( unsigned j = 0 ; j < nrays ; j++ ) {
	ofsreadings << readings[i][j] << "  ";
	ofsarea << area[i][j] << " ";
	ofsalpha << alpha[i][j] << "  " ;
	ofsalphavalid << ( alpha_valid[i][j] ? 1 : 0 ) << " ";
	ofsp3d << p3d[i][j] << " " ;
	ofsvalid << ( valid[i][j] ? 1 : 0 ) << " " ;
      }
      ofsreadings << endl;
      ofsarea << endl;
      ofsalpha << endl;
      ofsalphavalid << endl;
      ofsp3d << endl;
      ofsvalid << endl;
    }

    ofsreadings.close();
    ofsarea.close();
    ofsalpha.close();
    ofsalphavalid.close();
    ofsp3d.close();
    ofsvalid.close();
    

  }


}
