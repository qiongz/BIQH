#include"init.h"
void usage(char *target) {
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -n                       nphi\n";
    std::cout<<"  -e                       No. of electrons\n";
    std::cout<<"  -j                       J in [0,nphi-1]\n";
    std::cout<<"  -g                       aspect ratio gamma=lx/ly\n";
    std::cout<<"  -m                       Lambda\n";
    std::cout<<"Default: (nphi,n,lambda) = (2,2,200)\n";
}

void init_argv(int &nphi, int& nel,int &J, double & gamma,int &lambda,int argc,char *argv[])
{
    extern char *optarg;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"n:e:j:g:m:h:"))!=-1) {
        switch(ch) {
        case 'n':
            nphi=atoi(optarg);
            break;
        case 'e':
            nel=atoi(optarg);
            break;
        case 'j':
            J=atoi(optarg);
            break;
        case 'g':
            gamma=atof(optarg);
            break;
        case 'm':
            lambda=atoi(optarg);
            break;
        case 'h':
            errFlag++;
            break;
        default:
            errFlag++;
            break;
        }
    }
    try {
        if(nel<=0)
            throw std::logic_error("-n: positive value required !");
    } catch(std::logic_error &e) {
        std::cout<<e.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }
    catch(std::overflow_error &e) {
        std::cout<<e.what()<<std::endl;
        exit(2);
    }

    if(errFlag) {
        usage(argv[0]);
        exit(0);
    }
}
