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
#include<gsl/gsl_integration.h>
#include<gsl/gsl_sf.h>

#ifdef mkl
#define MKL_Complex16 std::complex<double>
#include"mkl.h"
#else
extern "C" int dsyevd_(char *, char *, int *, double *, int*, double *, double *, int *,int*,int*, int *);
extern "C" int zheevd_(char*, char*, int*,std::complex<double>*,int*,double*,std::complex<double>*,int*,double*,int*,int*,int*,int*);
#endif
#if __cplusplus > 199711L
#include<random>
#else
#include"mt19937-64.h"
#endif

using namespace std;
void diag_dsyevd(double *h, double *e, int l);
void diag_zheevd(complex<double> *h, double *e, int l);
double func_ExpInt(double t, void *params);
double Integrate_ExpInt(double z) ;
double func_BesselInt(double k, void *params);
double Integrate_BesselInt(double r,double d) ;

class Vec {
public:
    std::vector< complex<double> > value;
    long size;

    Vec();
    Vec(long _size);
    Vec(long _size,const complex<double> _init);
    Vec(const Vec & rhs);
    ~Vec();

    void assign(long _size, const complex<double> _init);
    void init_random(unsigned);
    void init_random(long,unsigned);
    void clear();
    double normed();
    void normalize();

    // operator overloading
    Vec & operator=(const Vec & rhs);
    Vec & operator-=(const Vec & rhs);
    Vec & operator+=(const Vec & rhs);
    Vec & operator*=(const double & rhs);
    Vec & operator*=(const complex<double> & rhs);
    Vec & operator/=(const double & rhs);
    Vec operator+(const Vec &);
    Vec operator-(const Vec &);
    Vec operator*(const double &);
    Vec operator*(const complex<double> &);
    Vec operator/(const double &);
    complex<double> operator*(const Vec &);
    friend ostream & operator<<(ostream & os, const Vec &);
};

class Mat {
public:
    // double compressed Sparse Row (CSR) Data Structure
    // outer_size for multi-threads storing, for which inner_indices are stored randomly
    std::vector<long> outer_starts,outer_size,inner_indices;
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
