#ifndef LANCZOS_HAMILTONIAN_H
#define LANCZOS_HAMILTONIAN_H
#include"matrix.h"
#include"basis.h"
#include<thread>
#include<mutex>
#include<gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
void swap(Vec *a,Vec *b,Vec *c);


class lhamil {
public:
    unsigned seed;  //!< Seed for RNGs
    long nHilbert;  //!< Hilbert space size
    long lambda;    //!< Lanczos update steps
    long nLL,nphi;
    double lx,ly,d;
    double E0,Ec;      //!< Ground state eigen energy
    // index: alpha*nphi*off_head*nphi*nphi+q_y*off_head*nphi*nphi+q_x*nphi*nphi+n*nphi+m
    vector<double> Coulomb_matrix; //!< store the Coulomb interaction matrix elements
    vector<complex<double> > FT;
    basis sector;
    Mat H;  //!< Hamiltonian matrix in CSR format
    //Mat O;  //!< Operator matrix in CSR format
    std::vector<double> norm; //!< Normalization coefficients vector in Lanczos update
    std::vector<double> overlap; //!< Overlap coefficients vector in Lanczos update
    std::vector< complex<double> > psir_0; //!< Ground state wave function in real-space
    std::vector< complex<double> > psir_1; //!< first excited wave function in real-space
    std::vector<double> psi_0; //!<Ground state eigenvector in Krylov subspace
    std::vector<double> psi_1; //!<first excited state eigenvector in Krylov subspace
    std::vector<double> psi_n0; //!<First element of eigenvectors in Krylov subspace
    std::vector<double> eigenvalues; //!< Eigenvalues
    lhamil();  //!< Empty constructor
    /**
     \param _H hamiltonian matrix
     \param _nHilbert Hilbert space size
     \param _lambda Lanczos update steps
     \param _seed Seed for RNGs
    */
    lhamil(const Mat & _H,long _nHilbert,long _lambda,unsigned _seed);  //!< Constructor with hamiltonian matrix as input
    /**
     \param _sector Basis sector
     \param _lambda Lanczos update steps
     \param _seed Seed for RNGs
    */
    lhamil(long _lambda,unsigned _seed); //!< Constructor with basis sector as input
    ~lhamil(); //!< Destructor
    const lhamil & operator=(const lhamil &);
    /** \param _sector Basis sector
    */
    void set_hamil(double _lx, double _ly, long _nphi, long _nLL,double _d,double _Delta_SAS,double _Delta_V,double _Delta_Z,double _theta_B,double theta_x,double theta_y,int nthread);  //!< Initialize hamiltonian matrix
    void peer_set_hamil(double,double,double,double,double,double,int,long,long);
    void Gram_Schmidt_orthogonalization(Vec &, int);
    void coeff_update(); //!< Lanczos update implemenation utilizing the Mat class
    void coeff_explicit_update(); //!< Lanczos update implemenation written in explicit arrays
    void coeff_update_wopt(vector< complex<double> > O_phi_0);
    void diag();  //!< Diagonalize the full Lanczos hamiltonian

    void eigenstates_reconstruction(); //!< Transform |psi_0> to |psir_0>
    double Coulomb_interaction(int alpha,int q_x, int q_y);
    double ground_state_energy();    //!< Ground state energy
    double first_excited_state_energy();    //!< 1st excited state energy
    double occupatation_number(int alpha,int j);
    double pseudospin_Sz();
    double Sz();
    double upper_Sz();
    double down_Sz();
    double pseudospin_Sx();
    double spinflip_tunneling();

    double spectral_function(double omega,double eta); //!< Spectral function with spin, continued fraction version
    void init_Coulomb_matrix(double);
    void print_hamil_CSR(); //!< print the hamiltonian matrix in the CSR format
    void print_hamil(int n); //!< print the full hamiltonian matrix
    void print_lhamil(int n);  //!< print the Lanczos hamiltonian matrix with first n x n elements
    void print_eigen(int n);  //!< print the first n eigenvalues
    void save_to_file(const char* filename);  //!< save object to file "filename"
    void read_from_file(const char*);        //!< load object from file "filename"
};
#endif
