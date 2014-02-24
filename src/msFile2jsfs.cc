/*-------------------------------------------------
 *  msFile2jsfs.cc
 *  Calculates the JSFS out of a ms output file.
 *  
 *  Author:     Paul R. Staab
 *  Email:      staab ( at ) bio.lmu.de
 *  Date:       2012-10-05
 *  Licence:    GNU GPLv3 (or later)
 *-------------------------------------------------
 */

#include <Rcpp.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

using namespace std;
static bool debug = false;

int L=0;
vector<int> MutOfSitePop1(0), MutOfSitePop2(0);
int individual = 0;
bool insideSite = 0;

void printjsfs(int *jsfs, const int &s1, const int &s2);
void parseLine(string line, std::vector<int> &jsfs, const int &s1, const int &s2);


RcppExport SEXP msFile2jsfs(SEXP filename_, 
                            SEXP s1_,
                            SEXP s2_) {

  BEGIN_RCPP    
  std::string filename = Rcpp::as<std::string>(filename_);
  int s1 = Rcpp::as<int>(s1_);
  int s2 = Rcpp::as<int>(s2_);

  std::vector<int> jsfs = std::vector<int>((s1+1)*(s2+1), 0);

  string line;
  ifstream myfile(filename.c_str());

  if (!myfile.is_open()) {
    throw std::runtime_error("Cannot open file");
  }

  while ( myfile.good() ) {
    getline(myfile, line);
    parseLine(line, jsfs, s1, s2);
  }   

  myfile.close();

  return Rcpp::wrap( jsfs );
  END_RCPP
}


void parseLine(string line, std::vector<int> &jsfs, const int &s1, const int &s2) {
  //if (debug) cout << line << endl;
  if (line == "") {
    //if (debug) cout << "> Empty line detected" << endl;
    if (insideSite){
      //Site ended
      insideSite = 0;

      //Calc JSFS
      for(int i=0; i<L; ++i){
        // if (debug) cout << MutOfSitePop1[i] << " : " << MutOfSitePop2[i] << endl;
        int index =  MutOfSitePop1[i] * (s2+1) + MutOfSitePop2[i];
        // if (debug) cout << "index:" << index << endl;
        if (index < 0 || index >= (s1+1)*(s2+1)) throw std::runtime_error("Wrong array index!");
        jsfs[ index ]++;
      }

      //printjsfs(jsfs);

      //Unset Site spezific variables
      individual = 0;
      L = 0;
      // if (debug) cout << "> Site ended detected" << endl;
    }
  }
  else if (line == "//") {
    insideSite = 1;
    //if (debug) cout << "> New Site detected" << endl;
  }
  else if ( (line.substr(0,1) == "0" || line.substr(0,1) == "1") && insideSite){
    //if (debug) cout << "> Data line detected" << endl;

    //Initialize per site variables on first data line
    if (L == 0) {
      L = line.length();
      MutOfSitePop1.resize(L);
      MutOfSitePop2.resize(L);
      for (int i=0; i<L; i++) {
        MutOfSitePop1[i] = 0;
        MutOfSitePop2[i] = 0;
      }
    }

    // if (debug) cout << "> Individual " << individual << "| Population " << individual / s1 << endl;

    //Parse the individual
    if (individual / s1 == 0){
      //we are in population 1 here
      for(int p=0; p<L; ++p) { 
        if (line[p] == '1') MutOfSitePop1[p]++;
      }

    }
    else {
      //we are in population 2 here
      for(int p=0; p<L; ++p) {
        if (line[p] == '1') MutOfSitePop2[p]++;
      }
    }

    individual++;
  }
}

/*void printjsfs(int *jsfs, const int &s1, const int &s2) {
  for (int i=0; i<=s1; i++) {
  for (int j=0; j<=s2; j++) {
//cout << jsfs[ i * (s2 + 1) + j ] << " ";
}
//cout << endl;
}
}*/

/*int main() {
  int s1=10, s2=15;
  int jsfs[176] = {0};
  int *pjsfs = jsfs;
  char filename[] = "msoutput";
  char *pfilename = filename;

  msFile2jsfs(&pfilename, &s1, &s2, jsfs);
  printjsfs(jsfs, s1, s2);

  return(0);
  }*/
