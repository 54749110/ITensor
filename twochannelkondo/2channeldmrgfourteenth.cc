#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include  "kondositeset.h"

using namespace itensor;
using channelkondo= twochannelkondosimpletwo<ElectronSite,SpinHalfSite>;

int 
main()
    {
// N must be 3Z
    int N = 11;
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

   auto t=1;
   auto Jk=0.1;
     for( int cycle= 0.1 ;cycle <= 120 ; cycle += 1 )  
    {
      Jk += 0.1 ;
    
    auto ampo = AutoMPO(sites);
     for (int j = 1; j <= (N-3)/2 ; j += 1) //electron hopping
     {
         ampo += -t, "Cdagup", j+1, "Cup", j;
         ampo += -t, "Cdagup", j,   "Cup", j+1;
         ampo += -t, "Cdagdn", j+1, "Cdn", j;
         ampo += -t, "Cdagdn", j,   "Cdn", j+1;
    }
    for(int j = 1; j <= N-1; j += 1) //2channel kondo hopping 1
    {
        ampo += 0.5 * Jk, "S+", N, "S-", j;
        ampo += 0.5 * Jk, "S-", N, "S+", j;
        ampo +=       Jk, "Sz", N, "Sz", j;
    }
        // ampo += 0.5 * Jk, "S+", N, "S-", (N-3)/4;
        // ampo += 0.5 * Jk, "S-", N, "S+", (N-3)/4;
        // ampo +=       Jk, "Sz", N, "Sz", (N-3)/4;
    //  for(int j = (N+3)/2; j <= N; j += 1) //2channel kondo hopping 1
    //  {
    //      ampo += 0.5 * Jk, "S+", (N+1)/2, "S-", j;
    //      ampo += 0.5 * Jk, "S-", (N+1)/2, "S+", j;
    //      ampo +=       Jk, "Sz", (N+1)/2, "Sz", j;
    //  }
    for (int j = (N+1)/2; j <= N-2 ; j += 1) //electron hopping
    {
        ampo += -t, "Cdagup", j+1, "Cup", j;
        ampo += -t, "Cdagup", j,   "Cup", j+1;
        ampo += -t, "Cdagdn", j+1, "Cdn", j;
        ampo += -t, "Cdagdn", j,   "Cdn", j+1;
    }
    // for (int j = 1; j <= N-2 ; j += 1) //electron hopping
    // {
    //     ampo += -t, "Cdagup", j+1, "Cup", j;
    //     ampo += -t, "Cdagup", j,   "Cup", j+1;
    //     ampo += -t, "Cdagdn", j+1, "Cdn", j;
    //     ampo += -t, "Cdagdn", j,   "Cdn", j+1;
    // }
    
    //     ampo += 0.5 * Jk, "S+", N, "S-", N-1;
    //     ampo += 0.5 * Jk, "S-", N, "S+", N-1;
    //     ampo +=       Jk, "Sz", N, "Sz", N-1;
    

    // PBC 1
        ampo += -t, "Cdagup",  (N-1)/2, "Cup", 1;
        ampo += -t, "Cdagup", 1,   "Cup", (N-1)/2;
        ampo += -t, "Cdagdn", (N-1)/2, "Cdn", 1;
        ampo += -t, "Cdagdn", 1,   "Cdn", (N-1)/2;
        ampo += -t, "Cdagup",  N-1, "Cup", (N+1)/2;
        ampo += -t, "Cdagup", (N+1)/2,   "Cup", N-1;
        ampo += -t, "Cdagdn", N-1, "Cdn", (N+1)/2;
        ampo += -t, "Cdagdn", (N+1)/2,   "Cdn", N-1;
    
    // // PBC 2
    //     ampo += -t, "Cdagup", N-1, "Cup", 1;
    //     ampo += -t, "Cdagup", 1,   "Cup", N-1;
    //     ampo += -t, "Cdagdn", N-1, "Cdn", 1;
    //     ampo += -t, "Cdagdn", 1,   "Cdn", N-1;

    H = toMPO(ampo);

    // Several tries to the occupation of initial state
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
    //print(s) ;
    


    // calculate average energy over kondo part
    // auto ampo_kondo = AutoMPO(sites);
    // for(int j = 1; j <= N-1; j += 1) //2channel kondo hopping 1
    // {
    //     ampo_kondo += 0.5 * Jk, "S+", N, "S-", j;
    //     ampo_kondo += 0.5 * Jk, "S-", N, "S+", j;
    //     ampo_kondo +=       Jk, "Sz", N, "Sz", j;
    // }
    // for(int j = 1; j <= (N-1)/2; j += 1) //2channel kondo hopping 1
    // {
    //     ampo_kondo += 0.5 , "S+", (N+1)/2, "S-", j;
    //     ampo_kondo += 0.5 , "S-", (N+1)/2, "S+", j;
    //     ampo_kondo +=  1  , "Sz", (N+1)/2, "Sz", j;
    // }
    // for(int j = (N+3)/2; j <= N; j += 1) //2channel kondo hopping 1
    // {
    //     ampo_kondo += 0.5 , "S+", (N+1)/2, "S-", j;
    //     ampo_kondo += 0.5 , "S-", (N+1)/2, "S+", j;
    //     ampo_kondo +=  1  , "Sz", (N+1)/2, "Sz", j;
    // }
    //  auto Hkondo = toMPO(ampo_kondo);
    //  printfln("\nKondo Energy = %.10f", inner(psi,Hkondo,psi) );
    

//for ( int s= 1 ;s <= 100 ; s += 1 ) 
//{
//print("\n", Groundstateenergy[s]) ;
//}
//    }
//Given an MPS called "psi",
//nd some particular bond "b" (1 <= b < length(psi))
//across which we want to compute the von Neumann entanglement


int b=N-1 ;
//"Gauge" the MPS to site b
psi.position(b);  

//SVD this wavefunction to get the spectrum
//of density-matrix eigenvalues
auto l = leftLinkIndex(psi,b);
auto s = siteIndex(psi,b);
auto [U,S,V] = svd(psi(b),{l,s});
auto u = commonIndex(U,S);

//Apply von Neumann formula
//to the squares of the singular values
Real SvN = 0.;
for(auto n : range1(dim(u)))
    {
    auto Sn = elt(S,n,n);
    auto p = sqr(Sn);
    if(p > 1E-12) SvN += -p*log(p);
    }
printfln("Across bond b=%d, SvN = %.10f",b,SvN);
print(Jk);
    }

    return 0;
    }