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
    /** Hamiltonian matrix in CSR format */
    std::vector< complex<double> > hamiltonian;
    /** Eigenvalues of the hamiltonian */
    std::vector<double> eigenvalues;
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
    void set_hamil(double _lx, double _ly, long _nphi, long _nLL,double _d,double _Delta_SAS,double _Delta_V,double theta_x, double theta_y,int nthread);  //!< Initialize hamiltonian matrix
    void peer_set_hamil(double,double,double,double,int,long,long);

    const hamil & operator=(const hamil &);
    /** Return the ground state energy of the system */
    double ground_state_energy();
    double pseudospin_Sz();
    double pseudospin_Sx();

    /** Diagonalize the full hamiltonian */
    void init_Coulomb_matrix(double theta_x);
    double Coulomb_interaction(int alpha,int q_x, int q_y);
    void diag();

    double spectral_function(vector< complex<double> >& O_phi_0,double omega,double _E0,double eta, int annil); //!< Spectral moments with spin
    /** Print the hamiltonian matrix */
    void print_hamil(int range);
    /** Print the eigenvalues of the system */
    void print_eigen(int range);
};
#endif
