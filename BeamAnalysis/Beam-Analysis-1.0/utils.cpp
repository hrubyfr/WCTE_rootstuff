
#include "utils.h"
#include <iostream>


int binary_search(const  double value, const std::vector< double> points){
    /*
        returns the index in which to place "value" in "points" to maintain order 

        if the object is provided 0.5 

        [0, 1, 2, 3]
    */

   if(value<points[0]){
    return -1;
   }else if(value>points[points.size()-1]){
    return points.size();
   }

   int min_abs = 0;
   int max_abs = points.size()-1;
   int lower_bin = abs(max_abs - min_abs)/2;
   lower_bin /= 2;
   int upper_bin = lower_bin+1; 

    while( !(points[lower_bin]<=value && points[upper_bin]>=value) ){
        if( abs(max_abs - min_abs)<=1){
            std::cout <<max_abs<<std::endl;
            throw std::runtime_error("Should be unreachable - are the points ordered?");
        }

        if(value<points[lower_bin]){
            max_abs = lower_bin;
        }
        if (value>points[upper_bin]){
            min_abs=upper_bin;
        }

        lower_bin = max_abs-min_abs;
        lower_bin/=2;
        lower_bin += min_abs;
        
        upper_bin = lower_bin+1;
    }

    if( !(value>=points[lower_bin] && value<=points[upper_bin]) ){
        throw std::runtime_error("Somehow the point is not between two other entries");
    }

    return lower_bin;    

}
