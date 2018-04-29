#ifndef INIT_H
#define INIT_H
#include<unistd.h>
#include<cstdlib>
#include<cmath>
#include<iostream>
#include<stdexcept>
void usage(char *);
void init_argv(int & nLL, int &nphi, int& nel, int &nel_up,int &J_up, int &J_down, int &kx_up, int &kx_down, double &d,double &gamma, int &lambda,int argc,char *argv[]);
class Timer
{
public:
    Timer() {
        clock_gettime(CLOCK_REALTIME, &beg_);
    }

    double elapsed() {
        clock_gettime(CLOCK_REALTIME, &end_);
        return (end_.tv_sec - beg_.tv_sec) +
               (end_.tv_nsec - beg_.tv_nsec)/1000000000.0;
    }

    unsigned long nanoseconds() {
        clock_gettime(CLOCK_REALTIME, &end_);
        return (end_.tv_nsec - beg_.tv_nsec);
    }

    void reset() {
        clock_gettime(CLOCK_REALTIME, &beg_);
    }

private:
    timespec beg_, end_;
};
#endif
