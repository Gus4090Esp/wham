//
//  histogram.h
//  wham
//
//  Created by Gustavo Espinoza on 6/28/22.
//

#ifndef histogram_h
#define histogram_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

using namespace std;

class Histogram{
public:
    int count;
    double start;
    double end;
    double width;
    int nbins;
    double *vals;
    double *histo;
    
    Histogram();
    
    Histogram(double start_val, double end_val, double binwidth);
    
    void cleanup();
    
    void add(double data);
    
    void output_normalized(std::fstream &file);
    
};

#endif /* histogram_h */
