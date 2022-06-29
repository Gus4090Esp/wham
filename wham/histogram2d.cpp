//
//  histogram2d.cpp
//  wham
//
//  Created by Gustavo Espinoza on 6/28/22.
//

#include "histogram2d.h"


Histogram2d::Histogram2d() {
   count = 0;
   xstart = -8;
   xend = 8;
   xwidth = 0.1;
   xnbins = ceil((xend - xstart)/xwidth);
   xvals = new double[xnbins];
   ystart = -8;
   yend = 8;
   ywidth = 0.1;
   ynbins = ceil((yend - ystart)/ywidth);
   yvals = new double[ynbins];
   histo = new double[xnbins*ynbins];
   for (int i = 0; i < xnbins; i++) {
      xvals[i] = xstart+(i+1)*xwidth/2.0;
   };
   for (int i = 0; i < ynbins; i++) {
      yvals[i] = ystart+(i+1)*ywidth/2.0;
   };
   for (int i = 0; i < xnbins*ynbins; i++) {
      histo[i] = 0.0;
   };
};

// create based on desired values
Histogram2d::Histogram2d(double xstart_val, double xend_val, double xbinwidth,
                         double ystart_val, double yend_val, double ybinwidth) {
   count = 0;
   xstart = xstart_val;
   xend = xend_val;
   xwidth = xbinwidth;
   xnbins = ceil((xend - xstart)/xwidth);
   xvals = new double[xnbins];
   ystart = ystart_val;
   yend = yend_val;
   ywidth = ybinwidth;
   ynbins = ceil((yend - ystart)/ywidth);
   yvals = new double[ynbins];
   histo = new double[xnbins*ynbins];
   for (int i = 0; i < xnbins; i++) {
      xvals[i] = xstart+(i+0.5)*xwidth;
   };
   for (int i = 0; i < ynbins; i++) {
      yvals[i] = ystart+(i+0.5)*ywidth;
   };
   for (int i = 0; i < xnbins*ynbins; i++) {
      histo[i] = 0.0;
   };
};


void Histogram2d::cleanup() {
   delete[] xvals;
   delete[] yvals;
   delete[] histo;
};

void Histogram2d::add(double xdata, double ydata) {
   int xbin, ybin;
   if (xdata > xstart && xdata < xend && ydata > ystart && ydata < yend) {
      xbin = floor((xdata-xstart)/xwidth);
      ybin = floor((ydata-ystart)/ywidth);
      histo[xbin*ynbins+ybin]++;
      count++;
   };
};

// output normalized histogram
void Histogram2d::output_normalized(std::fstream &file) {
   char line[1024];
   for (int i = 0; i < xnbins; i++) {
      for (int j = 0; j < ynbins; j++) {
         sprintf(line, "%8.4f %8.4f %12.8f\n", xvals[i], yvals[j], histo[i*ynbins+j]/double(count)/xwidth/ywidth);
         file << line;
      };
      file << endl;
   };
};






