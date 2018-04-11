#ifndef INIT_H
#define INIT_H
#include<unistd.h>
#include<cstdlib>
#include<cmath>
#include<iostream>
#include<stdexcept>
void usage(char *);
void init_argv(int &nphi, int& nel, int & J,double & gamma,int &lambda,int argc,char *argv[]);
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
