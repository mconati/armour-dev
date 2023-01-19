#ifndef HEADER_H
#define HEADER_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <fstream>
#include <cstring>
#include <boost/numeric/interval.hpp>
#include <cassert>
#include <vector>
#include <cstdint>

#define TYPE double
#define TYPE_MAX DBL_MAX
#define TYPE_MIN DBL_MIN

#define WARNING_PRINT printf
#define PRINT printf

namespace bn = boost::numeric;
namespace bi = bn::interval_lib;

using Interval = bn::interval<
        TYPE, 
        bi::policies<
            bi::save_state<bi::rounded_transc_std<TYPE> >,
            bi::checking_base<TYPE>
        > 
    >;


using std::vector;
using std::cout;
using std::endl;
using std::sort;
using std::swap;
using std::min;
using std::max;

#endif