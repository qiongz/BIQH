/************************************************
Basis generation library for Bilayer integer
quantum Hall systems
   @author qiongzhu
   19/03/2018
   last updated
   03/07/2018
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
    long nbasis;     // No. of basis for up/down-layer electrons
    int K,J;  // total sum of j for up/down-layers and kx for up/down-layers
    vector<int> popcount_table;
    vector<int> nel_up_table;
    std::unordered_map<unsigned long, long> basis_set; // basis set of up/down-layer electrons, I-J table
    vector<unsigned long> id;     // reversal table, J->I, Lin's Table is a 2D array
    vector<short> basis_C;
    explicit basis();
    basis(int _nphi,int _nel, int _nel_up);
    basis(int _nphi,int _nel, int _nel_up, int _J,int _K);
    const basis & operator=(const basis &);
    ~basis();
    unsigned long factorial(int,int);
    long common_divisor(int,int);
    void init();
    void init(int _nphi,int _nel,int _nel_up);
    void init(int _nphi,int _nel, int _nel_up,int _J,int _K);
    void generate(long,long,long,unsigned long);
    void generate_all_density(long,long,long,unsigned long);
    int get_sign(unsigned long c,int n,int m,int nt,int mt,int t);
    int get_sign(unsigned long c,int n, int nt);
    unsigned long translate(unsigned long c, int k, int &sign) ;
    unsigned long inv_translate(unsigned long c, int k, int &sign) ;
    void clear();
    void prlong();
};

#endif
