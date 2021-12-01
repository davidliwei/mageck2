/*
 *  RRA_main_entry.cpp
 *  Implementation of Robust Rank Aggregation (RRA)
 *
 *  Created by Han Xu and Wei Li.
 *  Modified by Wei Li.
 *  Copyright 2017 Dana Farber Cancer Institute. All rights reserved.
 *
 */
#include "math_api.h"
#include "words.h"
#include "rvgs.h"
#include "rngs.h"
#include "classdef.h"
#include "fileio.h"
#include "RRA.h"

//C++ functions
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
/*
 * Main entry
 */
int main (int argc, const char * argv[]) {

  vector<string> rra_argv;
  for(int i=0;i<argc;i++){
	  rra_argv.push_back(string(argv[i]));
  }
  
  return RRA_main(rra_argv);
  
}

