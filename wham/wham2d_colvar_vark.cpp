//
//  wham2d_colvar_vark.cpp
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

#include "histogram2d.h"

using namespace std;



void get_centers(std::ifstream &file, double centers1[], double centers2[], double ks1[], double ks2[], int nwindows);



int main(int narg, char **arg)
{
   if (narg != 2) {
      cout << "Need parameter input!!" << endl;
      exit(0);
   }

   double kt = 0.9;
   int nwindows = 32;
   char fname1[80];
   char coldir[80];
   char colnm[80];
   double binwidth1 = 0.02;
   double start_val1 = 24.02;
   double end_val1 = 31.8;
   double binwidth2 = 0.02;
   double start_val2 = 24.02;
   double end_val2 = 31.8;
   double tol = 0.00001;

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
      if (!strcmp(keyword, "num_windows")) {
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
      } else if (!strcmp(keyword, "bin_width2")) {
         binwidth2 = atof(value);
      } else if (!strcmp(keyword, "start_value2")) {
         start_val2 = atof(value);
      } else if (!strcmp(keyword, "end_value2")) {
         end_val2 = atof(value);
      } else if (!strcmp(keyword, "kbt")) {
         kt = atof(value);
      } else {
         cerr << "Could not recognize " << keyword << " keyword provided in "
              << arg[1] << endl;
      }
   }
   script.close();
   cout << "kt: " << kt << endl;

   double centers1[nwindows];
   double centers2[nwindows];
   double ks1[nwindows];  // window dependent spring constants
   double ks2[nwindows];  // window dependent spring constants
   f.open(fname1, fstream::in);
   get_centers(f, centers1, centers2, ks1, ks2, nwindows);
   f.close();

   Histogram2d hists[nwindows];

   cout << "Creating " << nwindows << " histograms" << endl;
   cout << "Using files named " << coldir << "/" << colnm << endl;

   double ddum, datax, datay;
   string sdum;
   int count = 0;
   fstream of;
   for (int i = 0; i < nwindows; i++) {
      count = 0;
      Histogram2d h(start_val1, end_val1, binwidth1, start_val2, end_val2, binwidth2);
      //cout << "reading window: " << i << endl;
      sprintf(fname1, "%s/%s", coldir, colnm);
      sprintf(fname1, fname1, i);
      f.open(fname1, fstream::in);
      while(!f.eof()) {
         f >> ddum >> datax >> datay;
         //cout << "line: " << count << endl;
         h.add(datax, datay);
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
   double Xvals[hists[0].xnbins];
   double Yvals[hists[0].ynbins];
   double NE[hists[0].xnbins*hists[0].ynbins];
   double Omega[hists[0].xnbins*hists[0].ynbins];
   double oldprobs[hists[0].xnbins*hists[0].ynbins];
   double maxdiff, diff;
   char line[1024];
   double min_p0 = 100.0;
   double probx[hists[0].xnbins];
   double proby[hists[0].ynbins];
   for (int i = 0; i < hists[0].xnbins; i++) {
      Xvals[i] = hists[0].xvals[i];
      for (int j = 0; j < hists[0].ynbins; j++) {
         NE[i*hists[0].ynbins+j] = 0.0;
      }
   }
   for (int i = 0; i < hists[0].ynbins; i++) {
      Yvals[i] = hists[0].yvals[i];
   }

   int ndx;
   for (int i = 0; i < nwindows; i++) {
      for (int j = 0; j < hists[0].xnbins; j++) {
         for (int k = 0; k < hists[0].ynbins; k++) {
            ndx = j*hists[0].ynbins+k;
            NE[ndx] += hists[i].histo[ndx];
         }
      }
   }

   for (int iter = 0; iter < 5000; iter++) {
      for (int i = 0; i < hists[0].xnbins*hists[0].ynbins; i++) {
         Omega[i] = 0.0;
      }
      for (int i = 0; i < nwindows; i++) {
         for (int j = 0; j < hists[0].xnbins; j++) {
            for (int k = 0; k < hists[0].ynbins; k++) {
               ndx = j*hists[0].ynbins+k;
               Omega[ndx] += hists[i].count
                           * exp(-ks1[i]*0.5/kt*(Xvals[j] - centers1[i])
                           * (Xvals[j] - centers1[i]))
                           * exp(-ks2[i]*0.5/kt*(Yvals[k] - centers2[i])
                           * (Yvals[k] - centers2[i])) / Z[i];
            }
         }
      }

      for (int i = 0; i < nwindows; i++) {
         Z[i] = 0.0;
         for (int j = 0; j < hists[0].xnbins; j++) {
            for (int k = 0; k < hists[0].ynbins; k++) {
               ndx = j*hists[0].ynbins+k;
               Z[i] += NE[ndx] * exp(-ks1[i]*0.5/kt*(Xvals[j] - centers1[i])
                               * (Xvals[j] - centers1[i]))
                               * exp(-ks2[i]*0.5/kt*(Yvals[k] - centers2[i])
                               * (Yvals[k] - centers2[i])) / Omega[ndx];
            }
         }
      }
      if (iter % 100 == 0) {
         // output Z's for convergence
         cout << "iteration: " << iter << " ";
         for (int i = 0; i < nwindows; i++) {
            //cout << centers1[i] << " " << centers2[i] << " " << Z[i] << endl;
            cout << Z[i] << " ";
         }
         cout << endl;
         min_p0 = 100.0;
         for (int i = 0; i < hists[0].xnbins*hists[0].ynbins; i++) {
            if (-log(NE[i]/Omega[i]) < min_p0) min_p0 = -log(NE[i]/Omega[i]);
         }
         sprintf(fname1, "wham2d_neglnp0_iter%06d.txt", iter);
         of.open(fname1, fstream::out);
         for (int i = 0; i < hists[0].xnbins; i++) {
            for (int j = 0; j < hists[0].ynbins; j++) {
               ndx = i*hists[0].ynbins+j;
               if (NE[ndx] > 0) {
                  sprintf(line, "%8.4f %8.4f %12.6f\n", Xvals[i], Yvals[j], -log(NE[ndx]/Omega[ndx]) - min_p0);
                  of << line;
               } else {
                  sprintf(line, "%8.4f %8.4f %s\n", Xvals[i], Yvals[j], "inf");
                  of << line;
               }
            }
            of << endl;
         }
         of.close();
         // output 1D free energies
         min_p0 = 1000.0;
         for (int i = 0; i < hists[0].xnbins; i++) {
            probx[i] = 0.0;
            for (int j = 0; j < hists[0].ynbins; j++) {
               ndx = i*hists[0].ynbins+j;
               probx[i] += NE[ndx]/Omega[ndx]*hists[0].xwidth;
            }
            if (-log(probx[i]) < min_p0) min_p0 = -log(probx[i]);
         }
         sprintf(fname1, "wham1d_neglnp0x_iter%06d.txt", iter);
         of.open(fname1, fstream::out);
         for (int i = 0; i < hists[0].xnbins; i++) {
            sprintf(line, "%8.4f %12.6f\n", Xvals[i], -log(probx[i]) - min_p0);
            of << line;
         }
         of.close();
         min_p0 = 1000.0;
         for (int i = 0; i < hists[0].ynbins; i++) {
            proby[i] = 0.0;
            for (int j = 0; j < hists[0].xnbins; j++) {
               ndx = j*hists[0].ynbins+i;
               proby[i] += NE[ndx]/Omega[ndx]*hists[0].ywidth;
            }
            if (-log(proby[i]) < min_p0) min_p0 = -log(proby[i]);
         }
         sprintf(fname1, "wham1d_neglnp0y_iter%06d.txt", iter);
         of.open(fname1, fstream::out);
         for (int i = 0; i < hists[0].ynbins; i++) {
            sprintf(line, "%8.4f %12.6f\n", Yvals[i], -log(proby[i]) - min_p0);
            of << line;
         }
         of.close();
      }
      // initialize oldprobs on 1st iteration
      if (iter == 0) {
         for (int i = 0; i < hists[0].xnbins; i++) {
            for (int j = 0; j < hists[0].ynbins; j++) {
               ndx = i*hists[0].ynbins+j;
               oldprobs[ndx] = NE[ndx]/Omega[ndx];
            }
         }
         continue;
      }
      // check for convergence of probabilities
      maxdiff = 0.0;
      for (int i = 0; i < hists[0].xnbins; i++) {
         for (int j = 0; j < hists[0].ynbins; j++) {
            ndx = i*hists[0].ynbins+j;
            diff = abs(oldprobs[ndx] - NE[ndx]/Omega[ndx]);
            if (diff > maxdiff) maxdiff = diff;
            oldprobs[ndx] = NE[ndx]/Omega[ndx];
         }
      }
      if (maxdiff < tol) {
         cout << "iter " << iter << " maxdiff: " << maxdiff << " less than tol " << tol << endl;
         break;
      }
   }

   // Output NE / Omega, WHAM results
   min_p0 = 100.0;
   for (int i = 0; i < hists[0].xnbins*hists[0].ynbins; i++) {
      if (-log(NE[i]/Omega[i]) < min_p0) min_p0 = -log(NE[i]/Omega[i]);
   }
   sprintf(fname1, "wham2d_neglnp0.txt");
   of.open(fname1, fstream::out);
   for (int i = 0; i < hists[0].xnbins; i++) {
      for (int j = 0; j < hists[0].ynbins; j++) {
         ndx = i*hists[0].ynbins+j;
         if (NE[ndx] > 0) {
            sprintf(line, "%8.4f %8.4f %12.6f\n", Xvals[i], Yvals[j], -log(NE[ndx]/Omega[ndx]) - min_p0);
            of << line;
         } else {
            sprintf(line, "%8.4f %8.4f %s\n", Xvals[i], Yvals[j], "inf");
            of << line;
         }
      }
      of << endl;
   }
   of.close();

   // output 1D free energies
   min_p0 = 1000.0;
   for (int i = 0; i < hists[0].xnbins; i++) {
      probx[i] = 0.0;
      for (int j = 0; j < hists[0].ynbins; j++) {
         ndx = i*hists[0].ynbins+j;
         probx[i] += NE[ndx]/Omega[ndx]*hists[0].xwidth;
      }
      if (-log(probx[i]) < min_p0) min_p0 = -log(probx[i]);
   }
   sprintf(fname1, "wham1d_neglnp0x.txt");
   of.open(fname1, fstream::out);
   for (int i = 0; i < hists[0].xnbins; i++) {
      sprintf(line, "%8.4f %12.6f\n", Xvals[i], -log(probx[i]) - min_p0);
      of << line;
   }
   of.close();
   min_p0 = 1000.0;
   for (int i = 0; i < hists[0].ynbins; i++) {
      proby[i] = 0.0;
      for (int j = 0; j < hists[0].xnbins; j++) {
         ndx = j*hists[0].ynbins+i;
         proby[i] += NE[ndx]/Omega[ndx]*hists[0].ywidth;
      }
      if (-log(proby[i]) < min_p0) min_p0 = -log(proby[i]);
   }
   sprintf(fname1, "wham1d_neglnp0y.txt");
   of.open(fname1, fstream::out);
   for (int i = 0; i < hists[0].ynbins; i++) {
      sprintf(line, "%8.4f %12.6f\n", Yvals[i], -log(proby[i]) - min_p0);
      of << line;
   }
   of.close();

   // cleanup
   for (int i = 0; i < nwindows; i++) {
      hists[i].cleanup();
   }

   return 0;
}

void get_centers(std::ifstream &file, double centers1[], double centers2[], double ks1[], double ks2[], int nwindows)
{
   for (int i = 0; i < nwindows; i++) {
      file >> centers1[i] >> centers2[i] >> ks1[i] >> ks2[i];
   }
}






