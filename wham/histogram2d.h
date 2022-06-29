//
//  histogram2d.h
//  wham
//
//  Created by Gustavo Espinoza on 6/28/22.
//

#ifndef histogram2d_h
#define histogram2d_h

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class Histogram2d {
    
public:
    int count;
    double xstart;
    double xend;
    double xwidth;
    int xnbins;
    double *xvals;
    double ystart;
    double yend;
    double ywidth;
    int ynbins;
    double *yvals;
    double *histo;
    
    
    Histogram2d();

    Histogram2d(double xstart_val, double xend_val, double xbinwidth,
                double ystart_yal, double yend_val, double ybinwidth);

    void cleanup();


    void add(double xdata, double ydata);

    
    void output_normalized(std::fstream &file);

 };


#endif /* histogram2d_h */
