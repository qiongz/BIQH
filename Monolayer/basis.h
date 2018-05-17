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
#include<unordered_map>
#include<vector>
#include<algorithm>
using namespace std;

class basis {
public:
    long nphi,nel,C,J,K;  // N_phi, No.of electrons, k_x, k_y
    long nbasis;     // No. of basis for electrons
    std::unordered_map<unsigned long,long> basis_set; // basis set of electrons, I-J table

    vector<unsigned long> id;     // reversal table, J->I, Lin's Table is a 2D array
    vector<int> popcount_table;
    vector<short> basis_C;
    explicit basis();
    basis(long _nphi,long _nel);
    basis(long _nphi,long _nel, long _J, long _K);
    const basis & operator=(const basis &);
    ~basis();
    long onsite_potential(long,long);
    unsigned long factorial(long,long);
    long common_divisor(long,long);
    void init();
    void init(long,long);
    void init(long,long,long,long);
    void generate(long,long,long,unsigned long);
    unsigned long translate(unsigned long,int,int &);
    unsigned long inv_translate(unsigned long,int,int &);
    long creation(long,long);
    long annihilation(long,long);
    int get_sign(long i ,long n ,long m,long nt,long mt);
    void prlong();
    void clear();
    void for_sum(long &a, int count,int index,int range, int n);
    friend ostream & operator<<(ostream & os, const basis &);
};

#endif
