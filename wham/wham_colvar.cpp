//
//  wham_colvar.cpp
//  wham
//
//  Created by Gustavo Espinoza on 6/28/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <vector>
#include <numeric>

#include "histogram.h"

using namespace std;



void get_centers(std::ifstream &file, double centers1[], int nwindows);

int main(int narg, char **arg)
{
   if (narg != 2) {
      cout << "Need parameter input!!" << endl;
      exit(0);
   }

   double kbias1 = 50.0;
   double kt = 0.9;
   int nwindows = 32;
   char fname1[80];
   char coldir[80];
   char colnm[80];
   double binwidth1 = 0.02;
   double start_val1 = 24.02;
   double end_val1 = 31.8;
   double tol = 0.000001;

   // read input script with parameters
   ifstream script, f;
   script.open(arg[1], fstream::in);
   char keyword[80];
   char value[80];
   string keyword_string;
   string value_string;

   while(script >> keyword_string) {
      script >> value_string;
      sprintf(keyword, "%s", keyword_string.c_str());
      sprintf(value, "%s", value_string.c_str());
      if (!strcmp(keyword, "kbias1")) {
         kbias1 = atof(value);
      } else if (!strcmp(keyword, "num_windows")) {
         nwindows = atoi(value);
      } else if (!strcmp(keyword, "centers_fname")) {
         sprintf(fname1, "%s", value);
      } else if (!strcmp(keyword, "colvar_dir")) {
         sprintf(coldir, "%s", value);
      } else if (!strcmp(keyword, "colvar_fname")) {
         sprintf(colnm, "%s", value);
      } else if (!strcmp(keyword, "bin_width1")) {
         binwidth1 = atof(value);
      } else if (!strcmp(keyword, "start_value1")) {
         start_val1 = atof(value);
      } else if (!strcmp(keyword, "end_value1")) {
         end_val1 = atof(value);
      } else if (!strcmp(keyword, "kbt")) {
         kt = atof(value);
      } else {
         cerr << "Could not recognize " << keyword << " keyword provided in "
              << arg[1] << endl;
      }
   }
   script.close();
   cout << "kbias1: " << kbias1 << " kt: " << kt << endl;

   double centers1[nwindows];
   //double ks[nwindows];  // window dependent spring constants
   f.open(fname1, fstream::in);
   get_centers(f, centers1, nwindows);
   f.close();

   Histogram hists[nwindows];

   cout << "Creating " << nwindows << " histograms" << endl;
   cout << "Using files named " << coldir << "/" << colnm << endl;

   double ddum, datax;
   string sdum;
   int count = 0;
   fstream of;
   for (int i = 0; i < nwindows; i++) {
      count = 0;
      Histogram h(start_val1, end_val1, binwidth1);
      //cout << "reading window: " << i << endl;
      sprintf(fname1, "%s/%s", coldir, colnm);
      sprintf(fname1, fname1, i);
      f.open(fname1, fstream::in);
      while(!f.eof()) {
         f >> ddum >> datax;
         //cout << "line: " << count << endl;
         h.add(datax);
         count++;
      }
      f.close();
      sprintf(fname1, "histogram%03d.txt", i);
      of.open(fname1, fstream::out);
      h.output_normalized(of);
      of.close();
      hists[i] = h;
   }

   cout << "WHAM: solving self-consistent equation" << endl;

   double Z[nwindows];
   for (int i = 0; i < nwindows; i++) {
      Z[i] = 1.0;
   }
   double Xvals[hists[0].nbins];
   double NE[hists[0].nbins];
   double Omega[hists[0].nbins];
   double oldprobs[hists[0].nbins];
   double maxdiff, diff;
   char line[1024];
   double min_p0 = 100.0;
   for (int i = 0; i < hists[0].nbins; i++) {
      Xvals[i] = hists[0].vals[i];
      NE[i] = 0.0;
   }

   for (int i = 0; i < nwindows; i++) {
      for (int j = 0; j < hists[0].nbins; j++) {
         NE[j] += hists[i].histo[j];
      }
   }

   for (int iter = 0; iter < 100000; iter++) {
      for (int i = 0; i < hists[0].nbins; i++) {
         Omega[i] = 0.0;
      }
      for (int i = 0; i < nwindows; i++) {
         for (int j = 0; j < hists[0].nbins; j++) {
            Omega[j] += hists[i].count
                        * exp(-kbias1*0.5/kt*(Xvals[j] - centers1[i])
                        * (Xvals[j] - centers1[i])) / Z[i];
         }
      }

      for (int i = 0; i < nwindows; i++) {
         Z[i] = 0.0;
         for (int j = 0; j < hists[0].nbins; j++) {
            Z[i] += NE[j] * exp(-kbias1*0.5/kt*(Xvals[j] - centers1[i])
                          * (Xvals[j] - centers1[i])) / Omega[j];
         }
      }
      if (iter % 1000 == 0) {
         // output Z's for convergence
         cout << "iteration: " << iter << " ";
         for (int i = 0; i < nwindows; i++) {
            cout << Z[i] << " ";
         }
         cout << endl;
         min_p0 = 100.0;
         for (int i = 0; i < hists[0].nbins; i++) {
            if (-log(NE[i]/Omega[i]) < min_p0) min_p0 = -log(NE[i]/Omega[i]);
         }
         sprintf(fname1, "wham_neglnp0_iter%06d.txt", iter);
         of.open(fname1, fstream::out);
         for (int i = 0; i < hists[0].nbins; i++) {
            sprintf(line, "%8.4f %12.6f\n", Xvals[i], -log(NE[i]/Omega[i]) - min_p0);
            of << line;
         }
         of.close();
      }
      // initialize oldprobs on 1st iteration
      if (iter == 0) {
         for (int i = 0; i < hists[0].nbins; i++) {
            oldprobs[i] = NE[i]/Omega[i];
         }
         continue;
      }
      // check for convergence of probabilities
      maxdiff = 0.0;
      for (int i = 0; i < hists[0].nbins; i++) {
         diff = abs(oldprobs[i] - NE[i]/Omega[i]);
         if (diff > maxdiff) maxdiff = diff;
         oldprobs[i] = NE[i]/Omega[i];
      }
      if (maxdiff < tol) {
         cout << "iter " << iter << " maxdiff: " << maxdiff << " less than tol " << tol << endl;
         break;
      }
   }

   // Output NE / Omega, WHAM results
   min_p0 = 100.0;
   for (int i = 0; i < hists[0].nbins; i++) {
      if (-log(NE[i]/Omega[i]) < min_p0) min_p0 = -log(NE[i]/Omega[i]);
   }
   sprintf(fname1, "wham_neglnp0.txt");
   of.open(fname1, fstream::out);
   for (int i = 0; i < hists[0].nbins; i++) {
      sprintf(line, "%8.4f %12.6f\n", Xvals[i], -log(NE[i]/Omega[i]) - min_p0);
      of << line;
   }
   of.close();

   // cleanup
   for (int i = 0; i < nwindows; i++) {
      hists[i].cleanup();
   }

   return 0;
}

void get_centers(std::ifstream &file, double centers1[], int nwindows)
{
   for (int i = 0; i < nwindows; i++) {
      file >> centers1[i];
   }
}
