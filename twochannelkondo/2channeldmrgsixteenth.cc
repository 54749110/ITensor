#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include  "kondositeset.h"

using namespace itensor;
using channelkondo= twochannelkondosimple<ElectronSite,SpinHalfSite>;

int 
main()
    {
// N must be 3Z
    int N =11;
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
      for( int cycle= 1 ;cycle <= 120 ; cycle += 1 )  
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
    for(int j = 1; j <= (N-1)/2; j += 1) //2channel kondo hopping 1
    {
        ampo += 0.5 * Jk, "S+", (N+1)/2, "S-", j;
        ampo += 0.5 * Jk, "S-", (N+1)/2, "S+", j;
        ampo +=       Jk, "Sz", (N+1)/2, "Sz", j;
    }
        // ampo += 0.5 * Jk, "S+", N, "S-", (N-3)/4;
        // ampo += 0.5 * Jk, "S-", N, "S+", (N-3)/4;
        // ampo +=       Jk, "Sz", N, "Sz", (N-3)/4;
     for(int j = (N+3)/2; j <= N; j += 1) //2channel kondo hopping 1
     {
         ampo += 0.5 * Jk, "S+", (N+1)/2, "S-", j;
         ampo += 0.5 * Jk, "S-", (N+1)/2, "S+", j;
         ampo +=       Jk, "Sz", (N+1)/2, "Sz", j;
     }
    for (int j = (N+3)/2; j <= N-1 ; j += 1) //electron hopping
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
        ampo += -t, "Cdagup",  N, "Cup", (N+3)/2;
        ampo += -t, "Cdagup", (N+3)/2,   "Cup", N;
        ampo += -t, "Cdagdn", N, "Cdn", (N+3)/2;
        ampo += -t, "Cdagdn", (N+3)/2,   "Cdn", N;
    
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

auto b=(N+1)/2 ;

    // calculate average energy over kondo part
    auto ampo_kondo1 = AutoMPO(sites);
    for(int j = 1; j <= b-1; j += 1) //2channel kondo hopping 1
    {
        ampo_kondo1 += 0.5 * Jk, "S+", b, "S-", j;
        ampo_kondo1 += 0.5 * Jk, "S-", b, "S+", j;
        ampo_kondo1 +=       Jk, "Sz", b, "Sz", j;
    }
     auto Hkondoleft = toMPO(ampo_kondo1);
     printfln("\nLeft Energy = %.10f", inner(psi,Hkondoleft,psi) );

     auto ampo_kondo2 = AutoMPO(sites);
    for(int j = b+1; j <= N; j += 1) //2channel kondo hopping 1
    {
        ampo_kondo2 += 0.5 * Jk, "S+", b, "S-", j;
        ampo_kondo2 += 0.5 * Jk, "S-", b, "S+", j;
        ampo_kondo2 +=       Jk, "Sz", b, "Sz", j;
    }
     auto Hkondoright = toMPO(ampo_kondo2);
     printfln("\nRight Energy = %.10f", inner(psi,Hkondoright,psi) );



psi.position(b); 

auto rhoright = psi(b) ;
for(int k=b+1 ; k<=N ;++k)
{
 rhoright *= psi(k) ;
}
rhoright =  prime(rhoright, "Site")*dag(rhoright) ;


auto T= inds(rhoright) ;
// Print(T) ;

auto R=T[0] ;
auto M=T[b] ;
auto W = M ;
 M=R ;
 R=W ;
rhoright =replaceInds(rhoright,{R,M},{M,R}) ;

auto [U,V] =diagHermitian(rhoright) ;
auto u = commonIndex(V,U);
double Svv = 0.000000;
for(auto n : range1(dim(u)))
    {
    auto Sn = elt(V,n,n);
  //  Print(Sn) ;
    Svv += abs(Sn) ;
    }
Svv = Svv-1 ;
Print(Svv) ;


psi.position(1) ;
auto rholeft = psi(1) ;
for(int k=2 ; k<= b ;++k)
{
 rholeft *= psi(k) ;
}
rholeft =  prime(rholeft, "Site")*dag(rholeft) ;


T= inds(rholeft) ;
// Print(T) ;

R=T[b-1] ;
M=T[N] ;
W = M ;
M=R ;
R=W ;
rholeft =replaceInds(rholeft,{R,M},{M,R}) ;

auto [U2,V2] = diagHermitian(rholeft) ;
u = commonIndex(V2,U2);
double Svv2 = 0.000000;
for(auto n : range1(dim(u)))
    {
    auto Sn2 = elt(V2,n,n);
  //  Print(Sn) ;
    Svv2 += abs(Sn2) ;
    }
Svv2 = Svv2-1 ;
Print(Svv2) ;

Print(Jk) ;
}
    return 0;
}