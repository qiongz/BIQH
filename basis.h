/************************************************
Basis generation library for Bilayer integer
quantum Hall systems
   @author qiongzhu
   19/3/2018
   Email:qiongzhunow@gmail.com
*************************************************/
#ifndef BASIS_H
#define BASIS_H
#include<iostream>
#include<string.h>
#include<iomanip>
#include<bitset>
#include<cstdlib>
#include<cmath>
#include<map>
#include<vector>
#include<algorithm>
using namespace std;

class basis {
public:
    long nphi,nel_up,nel_down, K;  // N_phi, up/down-layer electrons, total sum of k
    map<long,long> basis_up,basis_down; // basis set of up/down-layer electrons, I-J table

    long nbasis_up,nbasis_down;     // No. of basis for up/down-layer electrons
    vector<long> id_up,id_down;     // reversal table, J->I, Lin's Table is a 2D array
    explicit basis();
    basis(long _nphi,long _nel_up, long _nel_down);
    basis(long _nphi,long _nel_up, long _nel_down, long _K);
    const basis & operator=(const basis &);
    ~basis();
    // Delta_SAS, interlayer hopping
    long interlayer_hopping(long,long,long);
    long onsite_potential(long,long,long);
    long factorial(long,long);
    void init();
    void init(long,long,long);
    void generate_up(long);
    void generate_down(long);
    long creation(long,long);
    long annihilation(long,long);
    int get_signu(long i ,long n ,long m,long nt,long mt);
    int get_signd(long j ,long n ,long m,long nt,long mt);
    int get_signud(long i,long j,long n,long m,long nt,long mt);
    void prlong();
    friend ostream & operator<<(ostream & os, const basis &);
};

#endif
