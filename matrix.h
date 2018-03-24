#ifndef MAdoubleRIX_H
#define MAdoubleRIX_H
#include<cmath>
#include<complex>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<omp.h>
#include"mkl.h"

#if __cplusplus > 199711L
#include<random>
#else
#include"mt19937-64.h"
#endif
using namespace std;
//extern "C" int dsyev_(char *, char *, int *, double *, int*, double *, double *, int *, int *);
void diag_dsyev(double *h, double *e, int l);

class Vec {
public:
    std::vector<double> value;
    long size;

    Vec();
    Vec(long _size);
    Vec(long _size,const double _init);
    Vec(const Vec & rhs);
    ~Vec();

    void assign(long _size, const double _init);
    void init_random(unsigned);
    void init_random(long,unsigned);
    void clear();
    double normalize();

    // operator overloading
    Vec & operator=(const Vec & rhs);
    Vec & operator-=(const Vec & rhs);
    Vec & operator+=(const Vec & rhs);
    Vec & operator*=(const double & rhs);
    Vec & operator/=(const double & rhs);
    Vec operator+(const Vec &);
    Vec operator-(const Vec &);
    Vec operator*(const double &);
    Vec operator/(const double &);
    double operator*(const Vec &);
    friend ostream & operator<<(ostream & os, const Vec &);
};

class Mat {
public:
    // doublehe compressed Sparse Row (CSR) Data Structure
    std::vector<long> outer_starts,inner_indices;
    std::vector<double> value;

    Mat();
    Mat(const Mat &rhs);
    ~Mat();
    Mat & operator=(const Mat & rhs);
    // the last const means the object is a constant
    Vec operator*(const Vec &)const;
    vector<double> operator*(const vector<double> &)const;
    void init(const vector<long> &,const vector<long> &,const vector<double> &);
    void clear();
    void print();
};

#endif
