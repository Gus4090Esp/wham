//
//  histogram.cpp
//  wham
//
//  Created by Gustavo Espinoza on 6/28/22.
//

#include "histogram.h"

Histogram::Histogram() {
   count = 0;
   start = 20;
   end = 35;
   width = 0.1;
   nbins = ceil((end - start)/width);
   vals = new double[nbins];
   histo = new double[nbins];
   for (int i = 0; i < nbins; i++) {
      vals[i] = start+(i+1)*width/2.0;
      histo[i] = 0.0;
   };
};


Histogram::Histogram(double start_val, double end_val, double binwidth) {
   count = 0;
   start = start_val;
   end = end_val;
   width = binwidth;
   nbins = ceil((end - start)/width);
   vals = new double[nbins];
   histo = new double[nbins];
   for (int i = 0; i < nbins; i++) {
      vals[i] = start+(i+0.5)*width;
      histo[i] = 0.0;
   };
};


void Histogram::cleanup() {
   delete[] vals;
   delete[] histo;
};


void Histogram::add(double data) {
   int bin;
   if (data > start && data < end) {
      bin = floor((data-start)/width);
      histo[bin]++;
      count++;
   };
};

void Histogram::output_normalized(std::fstream &file) {
   char line[1024];
   for (int i = 0; i < nbins; i++) {
      sprintf(line, "%8.4f %12.8f\n", vals[i], histo[i]/double(count)/width);
      file << line;
   };
};
