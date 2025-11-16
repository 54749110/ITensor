#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    //Number of sites

    auto N = 5;

    auto sites = SpinHalf(N);

    //Make a random MPS for testing
    auto state = InitState(sites,"Up");
    auto psi = randomMPS(state);
    print(psi(2)) ;
    //auto D = prime(psi(2),"Site") ;
    //print(D) ;
    auto ir = commonIndex(psi(2),psi(3),"Link"); //index to right of bra
    auto C = (prime(prime(psi(2),"Site"),ir)) ;
    print(C) ;

    return 0;
    }