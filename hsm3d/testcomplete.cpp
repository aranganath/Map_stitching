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
#include "hough3d.h"
#include "io.h"

#include <iostream>
#include <vector>
#include <sys/time.h>
#include <time.h>
#include <fstream>
#include <ANN/ANN.h>
#include <getopt.h>

using namespace std;
using namespace h3d;
int dim =3;
istream* dataIn =NULL;
int npts=0;
int iter=0;
string inputFile;
bool readPt(istream &in, ANNpoint &p, point &pt)
{
  double value;
  for (int i = 0; i < dim; i++) {
    if(!(in >> value)) return false;
    p[i] = value;
    if(i==0) pt.x = value;
    if(i==1) pt.y = value;
    if(i==2) pt.z = value;
  }
  return true;
}

void printPt(ostream &out, ANNpoint &p){

  out << "(" << p[0];
  for (int i = 1; i < dim; i++) {
    out << ", " << p[i];
  }
  out << ")\n";
}

void processArgs(int argc, char **argv){
  const char* const short_opts = "d:";
  static ifstream datastream;
  static struct option long_opts[] = {"data-file", required_argument, nullptr,'d'};
  while(true){
    const auto opt = getopt_long(argc, argv, short_opts,long_opts,nullptr);
    if(-1==opt){
      break;
    }
    switch(opt){
      case 'd':
        datastream.open(string(optarg),ios::in);
        if(!datastream){
          cerr<<"Cannot open File!"<<endl;
          exit(1);
        }
        cout<<"Reading from file "<<string(optarg)<<endl;
        dataIn = &datastream;
        break;
        
      default:
        cerr << "Usage:\n\n"
        << "  testcomplete [-d data]"
        << "  where:\n"
        << "    data     name of file containing data points\n"
        << " Results are sent to the standard output.\n"
        << "\n"
        << " To run use:\n"
        << "    ./testcomplete -d <data pts file>\n";
        break;

    }
  }
}

int main(int argc,char **argv) {

  processArgs(argc,argv);
  ANNpointArray dataPts;
  ANNpoint queryPt;
  ANNidxArray nnIdx;
  ANNdistArray dists;
  queryPt = annAllocPt(dim);
  vector<point> PointSet;
  vector<point> newPointset;
  point it;
  dataPts = annAllocPts(300, dim);

  int npts = 0;
  while(readPt(*dataIn, dataPts[npts],it)){
    PointSet.push_back(it);
      printPt(cout, dataPts[npts]);
      npts++;
  }

  
  h3d_params p;
  
  p.max_norm = 1.1 * 50;
  
  matrix33 identity;
  point zero_point(0,0,0);
  
  h3d_struct ht1(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);
  h3d_struct ht2(p.angular_cell_size_deg,p.rho_cell_size,-p.max_norm,p.max_norm);
  
  ht1.compute_cube_ht(PointSet);
  
  new_ranked_solution randomrot;

  randomrot.rot.rotation.data[0][0] = 0.75;
  randomrot.rot.rotation.data[0][1] = -0.433;
  randomrot.rot.rotation.data[0][2] = 0.5;
  randomrot.rot.rotation.data[1][0] = 0.64951;
  randomrot.rot.rotation.data[1][1] = 0.625;
  randomrot.rot.rotation.data[1][2] = -0.433;
  randomrot.rot.rotation.data[2][0] = -0.125;
  randomrot.rot.rotation.data[2][1] = 0.6495190;
  randomrot.rot.rotation.data[2][2] = 0.75;
  randomrot.tran.direction.x = 20;
  randomrot.tran.direction.y = 20;
  randomrot.tran.direction.z = 20;
  for(int i=0;i<PointSet.size();++i){

    it = randomrot.rot.rotation * PointSet[i] + randomrot.tran.direction;
    newPointset.push_back(it);

  }

  ht2.compute_cube_ht(newPointset);
  list<rotation_hypothesis> hyp = ht1.guess_rotation(ht2);  



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
  cout<<"Evaluating rotation "<<randomrot.rot.rotation<<endl;
  cout<<"Evaluation translation"<<randomrot.tran.direction<<endl;
  cout << "Evaluating "  << roto_tran.size() << endl;
  vector<new_ranked_solution> v_final = ht1.evaluate_solutions_hausdorff_pts(newPointset, roto_tran, dataPts, npts);

  matrix33 m = randomrot.rot.rotation*v_final[0].rot.rotation;
  cout<<"Bestfit distance:"<<v_final[0].distance<<endl;
  cout<<"Bestfit rotation:"<<v_final[0].rot.rotation<<endl;
  cout<<"Bestfit translation:"<<v_final[0].tran.direction<<endl;
  cout<<"Secondfit distance:"<<v_final[1].distance<<endl;
  cout<<"Secondfit rotation:"<<v_final[1].rot.rotation<<endl;
  cout<<"Secondfit translation:"<<v_final[1].tran.direction<<endl;
  cout<<"Thirdfit distance:"<<v_final[2].distance<<endl;
  cout<<"Thirdfit rotation:"<<v_final[2].rot.rotation<<endl;
  cout<<"Thirdfit translation:"<<v_final[2].tran.direction<<endl;
  
  cout<<"Identity:"<<m<<endl;

  ofstream file;
  ofstream file1;
  file.open("FinalSolution.txt");
  file1.open("Original.txt");
  for(int i=0;i<newPointset.size();++i){
    matrix33 R = v_final[0].rot.rotation;
    point newP = R*(newPointset[i]-v_final[0].tran.direction);
    file<<newP.x<<" "<<newP.y<<" "<<newP.z<<endl;
    file1<<PointSet[i].x<<" "<<PointSet[i].y<<" "<<PointSet[i].z<<endl;
  }

  file1.close();
  file.close();
}
