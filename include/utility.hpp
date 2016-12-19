#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <iomanip>
#include <math.h>

std::vector<std::string> slice(std::vector<std::string>& inp,
                               size_t start,
                               size_t stop);
size_t indexof(std::vector<std::string>& v, std::string& val);

std::string sfloat(double v, unsigned int places);
std::string bp_formatter(unsigned int bp);
std::vector<double> arange(double start, double stop, double step);

inline double round(double val, unsigned int places) {  return roundf(val*pow(10, places)) / pow(10,places); } 
struct ValueRun {
    size_t start;
    size_t stop;
    int value;
    bool operator==(const ValueRun b);
    int length(void);
};

std::vector<ValueRun> runs_gte(const std::vector<int>& v, int thresh);
std::vector<ValueRun> runs_gte_classic(std::vector<int>& sequence, int minval, int minlength);

#endif
