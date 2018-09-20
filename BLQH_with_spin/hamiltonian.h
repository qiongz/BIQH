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
    void set_hamil(double _lx, double _ly, long _nphi, long _nLL,double _d,double _Delta_SAS,double _Delta_V,double _Delta_Z,int nthread);  //!< Initialize hamiltonian matrix
    void peer_set_hamil(double,double,double,int,long,long);

    const hamil & operator=(const hamil &);
    /** Return the ground state energy of the system */
    double ground_state_energy();
    double pseudospin_Sz();
    double pseudospin_Sx();
    double spinflip_tunneling();
    double Sz();

    /** Diagonalize the full hamiltonian */
    void init_Coulomb_matrix();
    double Coulomb_interaction(int alpha,int q_x, int q_y);
    void diag();

    /** Print the hamiltonian matrix */
    void print_hamil(int range);
    /** Print the eigenvalues of the system */
    void print_eigen(int range);
};
#endif
