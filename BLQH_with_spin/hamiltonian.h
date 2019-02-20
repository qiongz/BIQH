#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include"basis.h"
#include"matrix.h"
#include<thread>
#include<mutex>
class hamil {
public:
    /** Size of the Hilbert space */
    long nHilbert;
    /** seed for the RNGs */
    unsigned seed;
    long nphi,nLL;
    double lx,ly,d;
    double Ec;
    basis sector;
    // index: alpha*nphi*off_head*nphi*nphi+q_y*off_head*nphi*nphi+q_x*nphi*nphi+n*nphi+m
    std::vector<double> Coulomb_matrix; //!< store the Coulomb interaction matrix elements
    vector<complex<double> > FT;
    complex<double>* hamiltonian;
    /** Eigenvalues of the hamiltonian */
    double* eigenvalues;
    /** Ground state wave function */
    std::vector< complex<double> > psi_0;
    /** first excited state wave function */
    std::vector< complex<double> > psi_1;
    /** First element of all wave functions */
    std::vector< complex<double> > psi_n0;

    hamil();
    ~hamil();
    /** Initialize the hamiltonian matrix elements.
     \param sector basis sector,
     \param t hopping strength,
     \param U onsite replusive interaction strength
     */
    void set_hamil(const double &_lx, const double &_ly,const  long &_nphi, const long &_nLL,const double &_d,const double &_Delta_SAS,const double &_Delta_V,const double &_Delta_Z,const int &nthread);  //!< Initialize hamiltonian matrix
    void peer_set_hamil(const double&,const double&,const double&,const int&,const long&,const long&);

    const hamil & operator=(const hamil &);
    /** Return the ground state energy of the system */
    double ground_state_energy() const;
    double pseudospin_Sz()const;
    double pseudospin_Sx()const;
    double spinflip_tunneling()const;
    double Sz()const;

    /** Diagonalize the full hamiltonian */
    void init_Coulomb_matrix();
    double Coulomb_interaction(const int &alpha,const int &q_x,const  int &q_y);
    void diag();

    /** Print the hamiltonian matrix */
    void print_hamil(int &range)const;
    /** Print the eigenvalues of the system */
    void print_eigen(int &range)const;
};
#endif
