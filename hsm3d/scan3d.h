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


#ifndef __SCAN3D
#define __SCAN3D

#include "coordinates.h"

/*!
  \file scan3D.h
  \brief Provides scan3D, i.e. the type needed to store a three dimensional range scan
*/

namespace h3d {

  class h3d_struct; // forward declaration

  /*!
    \brief Represents the data returned by a three dimensional scanner
    
    It is assumed that the range scanner is associated with a righthand system
    with x pointing forward, z upwards, and therefore y to the right. It is implied
    that the range scanner rotates around its y axis.
  */
  class scan3D {

    friend class h3d_struct;
  public:
  

    /*! \brief Creates an instance with the given number of scans and beams per scan.

      After you have initialized readings and valid you must call scan_convert_points before calling compute_cube_ht.
      If you want to use the faster compute_cube_ht_normal you also have to call scan_find_normals.
      @param s number of scans
      @param b number of beams per scan */
    scan3D(unsigned int s,unsigned int b);

    /*! \brief Creates a dummy instance that cannot be used. Useful only if later on it is loaded from a file.*/
    scan3D();

    /*! \brief Deallocates memory */
    ~scan3D();

    /*! \brief Converts raw range readings into points */
    void scan_convert_points(void);

    /*! \brief Initalizes alpha and slpha_valid */
    void scan_find_normals(void);

    /*! \brief Loads a scan from the given json file.

      This function also calls scan_convert_points and scan_find_normals, so that the scan is ready to be processed.
      @param fname name of the file to load */
    int load_scan_jsonfile(char* fname);

    /*! \brief  Saves the scan to a series of text files. Useful only for internal debugging */
    void save_scan_to_file(const char*);

    /*! \brief Loads a scan from a file. Useful only for benchmarking and debug. */
    int load_absolute_points_from_file(const char*);

    int load_points_from_file(const char*);
    
    /*! \brief Getter method returning the number of horizontal scans in the 3D scan */
    unsigned int get_nscans() const { return nscans;  }

    /*! \brief Getter method returning the number of beams per horizontal scan */
    unsigned int get_nbeams() const { return nrays; }

    /*! */
    REAL get_tilt_deg(unsigned int i) const { return tilt_deg[i]; }

    /*! */
    REAL get_thetas_deg(unsigned int i) const { return thetas_deg[i]; }

    /*! \brief The actual range data retuned by the scanner, organized as a table with nscans rows and nbeams columns
     
      Public to avoid calling setter/getter methods. */
    REAL **readings;
    /*! \brief Boolean value indicating which beams are valid or not (for example because they are out of range) 
      
      Public to avoid calling setter/getter methods. */
    bool **valid;

  private:
    void allocate_memory();
    /*! \brief Vector defining the different rotation values in degrees */
    REAL *tilt_deg;
    /*! \brief Vector defining the orientation of the beams in the scannin plane (in degrees) */
    REAL *thetas_deg;
    /*! \brief Conversion of thetas_deg into radians (cached for efficiency reasons) */
    REAL *thetas_rad;
    /*! \brief Conversion of tilt_deg into radians (cached for efficiency reasons) */
    REAL *tilt_rad;
    /*! \brief Number of orizontal scans, i.e. number of elements in tilt_deg (or tilt_rad) */
    unsigned int nscans;
    /*! \brief Number of beams in a single scan, i.e. number of elements in thetas_deg (ord thetas_rad) */
    unsigned int nrays;
    /*!  \brief Location of the scanned points in the reference frame placed on the scanner (organized in a linear table). */
    point *points;
    /*! \brief Location of the scanned points in the reference fram placed on the scanner (kept in the rectangular structure). */
    point **p3d;
    /*! \brief Area surrounding a given point */
    REAL **area;
    /*! \brief Boolean flag indicating wether the normal is valid*/
    bool **alpha_valid;
    /*! \brief Ortogonal values to the planes defined by the beams */
    s2_point **alpha;
  };

}

#endif
