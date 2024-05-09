#pragma once

#include <iostream>
#include <cstdlib>
#include <signal.h>

// #define ASSERT(x) 

#define ASSERT(x) if (!(x)) { \
    std::cerr << "ASSERTION FAILED: [[ " << #x << " ]] at " << __FILE__ << ":" << __LINE__ << std::endl;\
    raise(SIGTRAP); \
}
