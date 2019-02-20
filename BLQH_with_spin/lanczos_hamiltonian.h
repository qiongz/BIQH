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
    lhamil(const Mat & _H,const long &_nHilbert,const long &_lambda,const unsigned &_seed);  //!< Constructor with hamiltonian matrix as input
    /**
     \param _sector Basis sector
     \param _lambda Lanczos update steps
     \param _seed Seed for RNGs
    */
    lhamil(const long &_lambda,const unsigned &_seed); //!< Constructor with basis sector as input
    ~lhamil(); //!< Destructor
    const lhamil & operator=(const lhamil &);
    /** \param _sector Basis sector
    */
    void set_hamil(const double &_lx, const double &_ly, const long &_nphi,const long &_nLL,const double &_d,const double &_Delta_SAS, const double &_Delta_V,const double &_Delta_Z,const double &_theta_B,const double &theta_x,const double &theta_y,const int &nthread);  //!< Initialize hamiltonian matrix
    void peer_set_hamil(const double&,const double&,const double&,const double&,const double &,const double&,const int&,const long&,const long&);
    void Gram_Schmidt_orthogonalization(Vec &, int);
    void coeff_update(); //!< Lanczos update implemenation utilizing the Mat class
    void coeff_explicit_update(); //!< Lanczos update implemenation written in explicit arrays
    void coeff_update_wopt(vector< complex<double> > O_phi_0);
    void diag();  //!< Diagonalize the full Lanczos hamiltonian

    void eigenstates_reconstruction(); //!< Transform |psi_0> to |psir_0>
    double Coulomb_interaction(const int &alpha,const int &q_x, const int &q_y);
    double ground_state_energy()const;    //!< Ground state energy
    double first_excited_state_energy()const;    //!< 1st excited state energy
    double occupatation_number(const int &alpha,const int &j)const;
    double pseudospin_Sz()const;
    double Sz()const;
    double upper_Sz()const;
    double down_Sz()const;
    double pseudospin_Sx()const;
    double spinflip_tunneling()const;

    double spectral_function(const double &omega,const double &eta); //!< Spectral function with spin, continued fraction version
    void init_Coulomb_matrix(const double&);
    void print_hamil_CSR()const; //!< print the hamiltonian matrix in the CSR format
    void print_hamil(int n)const; //!< print the full hamiltonian matrix
    void print_lhamil(int n)const;  //!< print the Lanczos hamiltonian matrix with first n x n elements
    void print_eigen(int n)const;  //!< print the first n eigenvalues
    void save_to_file(const char* filename);  //!< save object to file "filename"
    void read_from_file(const char*);        //!< load object from file "filename"
};
#endif
