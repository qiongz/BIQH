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
#include<unordered_map>
#include<climits>
#include<vector>
#include<algorithm>
using namespace std;

class basis {
public:
    int nphi,nel,nel_up,nel_down;  //N_phi, up/down-layer electrons
    int K,J,C;  // total sum of j for up/down-layers and kx for up/down-layers
    vector<int> popcount_table;
    std::unordered_map<unsigned long long, long> basis_set; // basis set of up/down-layer electrons, I-J table

    long nbasis;     // No. of basis for up/down-layer electrons
    vector<unsigned long long> id;     // reversal table, J->I, Lin's Table is a 2D array
    explicit basis();
    basis(int _nphi,int _nel_up, int _nel_down);
    basis(int _nphi,int _nel_up, int _nel_down,int _J,int _K);
    const basis & operator=(const basis &);
    ~basis();
    // Delta_SAS, interlayer hopping
    unsigned long factorial(int,int);
    long common_divisor(int,int);
    void init();
    void clear();
    void init(int _nphi,int _nel_up,int _nel_down);
    void init(int _nphi,int _nel_up, int _nel_down,int _J,int _K);
    void generate(long,long,long,unsigned long long);
    unsigned long long translate(unsigned long long,int,int &);
    unsigned long long inv_translate(unsigned long long,int,int &);
    //long creation(long,long);
    //long annihilation(long,long);
    int get_sign(unsigned long long c ,int n ,int m,int nt,int mt);
    void prlong();
    friend ostream & operator<<(ostream & os, const basis &);
};

#endif
