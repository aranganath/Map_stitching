#include <iostream>
#include <stdlib.h>
#include <fstream>
using namespace std;
int main(){

	ofstream myfile;
	myfile.open("PointTest.txt");
	for(int i=0; i<170;++i){

		float x = rand()%5;
		float y = rand()%5;
		float z = rand()%5;
		// cout<<x<<" "<<y<<" "<<z<<endl; 
		myfile<<x<<" "<<y<<" "<<z<<endl;

	}
	myfile.close();

	return 0;
}
