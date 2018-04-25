/************************************************
Basis generation library for Bilayer integer
quantum Hall systems
   @author qiongzhu
   19/3/2018
   Email:qiongzhunow@gmail.com
*************************************************/
#ifndef BASIS_H
#define BASIS_H
#include<cmath>
#include<iostream>
#include<string.h>
#include<iomanip>
#include<bitset>
#include<cstdlib>
#include<map>
#include<vector>
#include<algorithm>
using namespace std;

class basis {
public:
    long nphi,nel,C,J,K;  // N_phi, No.of electrons, k_x, k_y
    map<long,long> basis_set; // basis set of electrons, I-J table

    long nbasis;     // No. of basis for electrons
    vector<long> id;     // reversal table, J->I, Lin's Table is a 2D array
    explicit basis();
    basis(long _nphi,long _nel);
    basis(long _nphi,long _nel, long _J, long _K);
    const basis & operator=(const basis &);
    ~basis();
    long onsite_potential(long,long);
    long factorial(long,long);
    long common_divisor(long,long);
    void init();
    void init(long,long);
    void init(long,long,long,long);
    void generate(long,long,long,long);
    long translate(long,long);
    long inv_translate(long,long);
    long creation(long,long);
    long annihilation(long,long);
    int get_sign(long i ,long n ,long m,long nt,long mt);
    void prlong();
    void clear();
    void for_sum(long &a, int count,int index,int range, int n);
    friend ostream & operator<<(ostream & os, const basis &);
};

#endif
