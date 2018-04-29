#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include"basis.h"
#include"matrix.h"
class hamil {
public:
    /** Size of the Hilbert space */
    long nHilbert;
    /** seed for the RNGs */
    unsigned seed;
    long nphi,nLL;
    double lx,ly,d;
    double E_cl;
    // index: alpha*nphi*off_head*nphi*nphi+q_y*off_head*nphi*nphi+q_x*nphi*nphi+n*nphi+m
    vector<double> Coulomb_matrix; //!< store the Coulomb interaction matrix elements
    /** Hamiltonian matrix in CSR format */
    vector< complex<double> > hamiltonian;
    /** Eigenvalues of the hamiltonian */
    std::vector<double> eigenvalues;
    /** Ground state wave function */
    std::vector< complex<double> > psi_0;
    /** First element of all wave functions */
    std::vector< complex<double> > psi_n0;

    hamil();
    ~hamil();
    /** Initialize the hamiltonian matrix elements.
     \param sector basis sector,
     \param t hopping strength,
     \param U onsite replusive interaction strength
     */
    void set_hamil(basis & sector, double _lx, double _ly, long nphi,long nLL, double _d);
    const hamil & operator=(const hamil &);
    /** Return the ground state energy of the system */
    double ground_state_energy();
    /** Diagonalize the full hamiltonian */
    void init_Coulomb_matrix();
    double Coulomb_interaction(int alpha,int q_x, int q_y);
    void diag();

    double spectral_function(vector< complex<double> >& O_phi_0,double omega,double _E0,double eta, int annil); //!< Spectral moments with spin
    /** Print the hamiltonian matrix */
    void print_hamil(int range);
    /** Print the eigenvalues of the system */
    void print_eigen(int range);
};
#endif
