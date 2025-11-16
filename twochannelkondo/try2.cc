#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include  "kondositeset.h"

using namespace itensor;
using channelkondo= twochannelkondosimple<ElectronSite,SpinHalfSite>;

int 
main()
    {
// N must be 3Z
    int N = 3;
    auto j = 3 ;
    channelkondo sites;
    MPO H;
    MPS psi0;
    //
    // Initialize the site degrees of freedom
    // Setting "ConserveQNs=",true makes the indices
    // carry Sz quantum numbers and will lead to 
    // block-sparse MPO and MPS tensors
    //
    sites = channelkondo(N); //make a chain of N spin 1/2's

auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    psi0 = MPS(state);
    print(psi0) ;
    psi0.position(3) ;
    auto S= siteIndex(psi0,j) ;
    auto M= leftLinkIndex(psi0, j) ;
    printfln("INDS",S) ;
}