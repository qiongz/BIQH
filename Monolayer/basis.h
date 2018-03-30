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
    long nphi,nel, K;  // N_phi, No.of electrons, total sum of k
    map<long,long> basis_set; // basis set of electrons, I-J table

    long nbasis;     // No. of basis for electrons
    vector<long> id;     // reversal table, J->I, Lin's Table is a 2D array
    explicit basis();
    basis(long _nphi,long _nel);
    basis(long _nphi,long _nel, long _K);
    const basis & operator=(const basis &);
    ~basis();
    long onsite_potential(long,long);
    long factorial(long,long);
    void init();
    void init(long,long);
    void generate(long);
    long creation(long,long);
    long annihilation(long,long);
    int get_sign(long,long,long,long,long);
    void prlong();
    friend ostream & operator<<(ostream & os, const basis &);
};

#endif
