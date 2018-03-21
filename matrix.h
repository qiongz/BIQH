#ifndef MATRIX_H
#define MATRIX_H
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

#if __cplusplus > 199711L
#include<random>
#else
#include"mt19937-64.h"
#endif
using namespace std;

extern "C" int zheev_(char *, char *, int *,complex<double>*,int *, double *,complex<double>*, int *, double *, int *);
void diag_zheev(complex<double> *hamiltonian, double *energy, int l);

class Vec {
public:
    std::vector< complex<double> > value;
    long size;

    Vec();
    Vec(long _size);
    Vec(long _size,const double _init);
    Vec(const Vec & rhs);
    ~Vec();

    void assign(long _size, const complex<double> _init);
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
    std::vector< complex<double> > value;

    Mat();
    Mat(const Mat &rhs);
    ~Mat();
    Mat & operator=(const Mat & rhs);
    // the last const means the object is a constant
    Vec operator*(const Vec &)const;
    vector< complex<double> > operator*(const vector< complex<double> > &)const;
    void init(const vector<long> &,const vector<long> &,const vector< complex<double> > &);
    void clear();
    void print();
};

#endif
