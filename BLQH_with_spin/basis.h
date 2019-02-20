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
    int nphi,nel,nel_up,nel_down,nel_su,nel_sd;  //N_phi, up/down-layer, spin-up/down electrons
    long nbasis;     // No. of basis for up/down-layer electrons
    int K,J;  // total sum of j for up/down-layers and kx for up/down-layers
    vector<int> popcount_table;
    std::unordered_map<unsigned long,long> basis_set; // basis set of up/down-layer electrons, I-J table
    vector<unsigned long> id;     // reversal table, J->I, Lin's Table is a 2D array
    vector<short> basis_C;
    explicit basis();
    basis(const int & _nphi,const int & _nel);
    basis(const int & _nphi,const int & _nel, const int & _J,const int &_K);
    basis(const int & _nphi,const int & _nel, const int & _nel_up,const int & _J,const int &_K);
    const basis & operator=(const basis &);
    ~basis();
    void init();
    void init(const int &_nphi,const int &_nel);
    void init(const int &_nphi,const int & _nel,const int& _J,const int& _K);
    void init(const int &_nphi,const int &_nel,const int &_nel_up,const int& _J,const int &_K);
    void init(const int &_nphi,const int &_nel,const int &_nel_up,const int &_nel_su,const int &_J,const int &_K);
    void generate(long,long,long,unsigned long);
    void generate_all_density(long,long,long,unsigned long);
    void generate_fixed_density(long,long,long,unsigned long);
    int get_sign(const unsigned long & c,const int &n,const int &m,const int &nt,const int &mt,const int &t) const;
    int get_sign(const unsigned long & c,const int &n,const int &nt) const;
    int get_nel(const int &,const long &) const;
    int get_nel(const int &,const unsigned long &) const;
    unsigned long translate(const unsigned long &c, const int &k, int &sign) const;
    unsigned long inv_translate(const unsigned long &c, const int &k, int &sign) const;
    unsigned long factorial(const int &, const int &) const;
    long common_divisor(const int&,const int &) const ;
    void prlong() const;
    void clear();
};

#endif
