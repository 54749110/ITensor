#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include  "kondositeset.h"

using namespace itensor;
using channelkondo= twochannelkondo<SpinHalfSite,ElectronSite>;

int 
main()
    {
    int N = 10;
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

   //double Groundstateenergy[100] ;

   auto  t=1;
   auto Jk=0;
   for( int s= 1 ;s <= 10 ; s += 1 )  
    {
      Jk += 1 ;
    
    auto ampo = AutoMPO(sites);
    for (int j = 2; j <= N-2 ; j += 1) //electron hopping
    {
        ampo += -t, "Cdagup", j+1, "Cup", j;
        ampo += -t, "Cdagup", j,   "Cup", j+1;
        ampo += -t, "Cdagdn", j+1, "Cdn", j;
        ampo += -t, "Cdagdn", j,   "Cdn", j+1;
    }
    for(int j = 1; j <= N-2; j += 1) //2channel kondo hopping
    {
        ampo += 0.5 * Jk, "S+", 1, "S-", j+1;
        ampo += 0.5 * Jk, "S-", 1, "S+", j+1;
        ampo +=       Jk, "Sz", 1, "Sz", j+1;
        ampo += 0.5 * Jk, "S+", N, "S-", j+1;
        ampo += 0.5 * Jk, "S-", N, "S+", j+1;
        ampo +=       Jk, "Sz", N, "Sz", j+1;
    }
    H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    psi0 = MPS(state);

    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", inner(psi0,H,psi0) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
    //Groundstateenergy[s] =energy ;
    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing inner = %.10f", inner(psi,H,psi) );
    print(s) ;
    }


    //for ( int s= 1 ;s <= 100 ; s += 1 ) 
    //{
    //print("\n", Groundstateenergy[s]) ;
    //}

    return 0;
    }