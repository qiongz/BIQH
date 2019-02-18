#ifndef INIT_H
#define INIT_H
#include<unistd.h>
#include<cstdlib>
#include<cmath>
#include<iostream>
#include<stdexcept>
#if __cplusplus > 199711L
#include<chrono>
#endif
void usage(char *);
void init_argv(int & nLL, int &nphi, int& nel,int &nel_up,int &nel_su, int &J, int &kx, double &d,double & Delta_SAS,double &Delta_V,double &Delta_Z, double &gamma, int &lambda,double & theta,int &nthread,unsigned long &,int argc,char *argv[]);
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
