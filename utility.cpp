#include "utility.hpp"


// Common string and indexing functions


std::vector<std::string> slice(std::vector<std::string>& inp, size_t start, size_t stop) {
    std::vector<std::string> out;
    for (size_t i = start; i < stop; ++i) {
        out.push_back(inp[i]);
    }

    return out;
}



size_t indexof(std::vector<std::string>& v, std::string& val) {
    std::string tok;
    for (size_t i = 0; i < v.size(); ++i) {
        tok = v[i];
        if (!tok.compare(val)) {
            return i;
        }
    }
    return -1;
}

std::string sfloat(double v, unsigned int places) {
    std::stringstream s;
    s << std::fixed << std::setprecision(places) << v;
    return s.str(); 
}

std::string bp_formatter(unsigned int bp) { 
    if (bp > 1e9) return sfloat(bp/1e9, 2) + "Gb";
    if (bp > 1e6) return sfloat(bp/1e6, 2) + "Mb";
    if (bp > 1e3) return sfloat(bp/1e3, 2) + "kb";
    return std::to_string(bp);
}

std::vector<double> arange(double start, double stop, double step) {
    std::vector<double> out;
    
    double v = start; 
    while (v <= stop) {
        out.push_back(v);
        v += step;
    }

    return out;
}


int ValueRun::length(void) {
    return stop - start;
}
bool ValueRun::operator==(const ValueRun b) {
    return ((start == b.start) && (stop == b.stop) && (value == b.value));
}

std::vector<ValueRun> runs_gte(const std::vector<int>& v, int thresh) {
    std::vector<ValueRun> outp;
    size_t start, stop;
    bool in_run = false;
    for (size_t i = 0; i < v.size(); ++i) {
        int val = v[i];
        if (!in_run && (val >= thresh)) {
            start = i;
            in_run = true;
        } else if (in_run && (val < thresh)) {
            stop = i;
            ValueRun run = {start, stop, v[i-1]};
            outp.push_back(run);
            in_run = false;
        }
    }


    return outp;
}

std::vector<ValueRun> runs_gte_classic(std::vector<int>& sequence, int minval, int minlength) {
    std::vector<ValueRun> out;
    bool inrun = false;
    size_t start, stop, i;
    int v;

    i = 0;
    while (i < sequence.size()) {
        v = sequence[i];

        if ((!inrun) && (v >= minval)) {
            inrun = true;
            start = i;
        } else if (inrun && (v < minval)) {
            inrun = false;
            stop = i - 1;
            if ((stop - start) >= minlength) {
                out.push_back(ValueRun{start,stop,sequence[stop-1]});
            } 
        }

        i++;
    }
    if (inrun && (i-start) > minlength) {
        out.push_back(ValueRun{start, i, sequence[i-1]});
    }

    return out;
}

std::string current_time_string(void) {
    time_t rt = time(NULL); 
    struct tm * local;
    local = localtime(&rt);

    std::string s = asctime(local);
    
    return s;
}

std::string print_elapsed(const timeval& t) {
    const int seconds_per_minute = 60;
    const int seconds_per_hour = 60 * seconds_per_minute;
    const int seconds_per_day = seconds_per_hour * 24;
    
    long vals[4] = {0,0,0,0};
    const char* units = "dhms";    
    long seconds = t.tv_sec;
    long millisec = lround(t.tv_usec / 1000.0);

    ldiv_t quot = ldiv(seconds, seconds_per_day);
    vals[0] = quot.quot;
    seconds = quot.rem;
    
    quot = ldiv(seconds, seconds_per_hour);
    vals[1] = quot.quot;
    seconds = quot.rem;
    
    quot = ldiv(seconds, seconds_per_minute);
    vals[2] = quot.quot;
    vals[3] = quot.rem;



    int startunit = -1;

    for (int i = 0; i < 4; i++) {
      if (vals[i]) { 
        startunit = i;
        break;
      }
    }

    startunit = startunit >= 0 ? startunit : 3;


    std::stringstream ss; 
    for (int i = startunit; i < 3; i++) {
      ss <<  vals[i] << units[i];      
    }

    ss << vals[3] << '.' << millisec << units[3];
    return ss.str();
}

