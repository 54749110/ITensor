#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include  "kondositeset.h"

using namespace itensor;

int 
main()
    {
// N must be 6Z
// construct siteset
   int N = 3*96;
   MPO H;
   MPS psi0;
   SpinHalf sites ;
   sites = SpinHalf(N,{"ConserveSz",false}); //make a chain of N spin 1/2's


   auto J=4;
   auto h=5.1;
   auto g=0.5;

    
    // construct hamiltonian
    auto ampo = AutoMPO(sites);

     // chain 1 Ising OBC "Sz"=1/2*sigma_z
     for (int j = 1; j <= N-4 ; j += 3) 
     {
        ampo += -1*4, "Sz", j, "Sz", j+3 ;
     }
     for (int j = 1; j <= N-1 ; j += 3) 
     {
        ampo += -1, "S+", j ;
        ampo += -1, "S-", j ;
     }

    //chain 1 Ising PBC 
       ampo += -1*4, "Sz", N-2, "Sz", 1 ;

    //chain 2 Ising OBC
    for(int j = 2; j <= N; j += 3) 
    {
        ampo += 4 * J, "Sz", j, "Sz", j+1 ;
    }
    for(int j = 3; j <= N-1; j += 3) 
    {
        ampo += 4 * J, "Sz", j, "Sz", j+2 ;
    }
    for(int j = 2; j <= N; j += 3) 
    {
        ampo += -h, "S+" , j ;
        ampo += -h, "S-" , j ;
        ampo += -h, "S+" , j+1 ;
        ampo += -h, "S-" , j+1 ;
    }

    // chain 2 PBC
        ampo += 4 * J, "Sz", N, "Sz", 2 ;
    
    // Hopping 
    for (int j = 1; j <= N-1 ; j += 3) 
    {
       ampo += g*2, "S+", j, "Sz", j+1 ; 
       ampo += g*2, "S-", j, "Sz", j+1 ; 
    }
    for (int j = 1; j <= N-4 ; j += 3) 
    {
      ampo += -g*8, "Sz", j, "Sz", j+3 ,"Sz" , j+2 ;
    }
    
    // Hopping PBC
      ampo += -g*8, "Sz", N-2, "Sz", 1 ,"Sz" , N ;


    // // chain 1 Ising OBC "Sz"=1/2*sigma_z
    //  for (int j = 1; j <= N-4 ; j += 3) 
    //  {
    //     ampo += -1, "Sz", j, "Sz", j+3 ;
    //  }
    //  for (int j = 1; j <= N-1 ; j += 3) 
    //  {
    //     ampo += -0.5, "S+", j ;
    //     ampo += -0.5, "S-", j ;
    //  }

    // //chain 1 Ising PBC 
    //    ampo += -1, "Sz", N-2, "Sz", 1 ;

    // //chain 2 Ising OBC
    // for(int j = 2; j <= N; j += 3) 
    // {
    //     ampo +=  J, "Sz", j, "Sz", j+1 ;
    // }
    // for(int j = 3; j <= N-1; j += 3) 
    // {
    //     ampo +=  J, "Sz", j, "Sz", j+2 ;
    // }
    // for(int j = 2; j <= N; j += 3) 
    // {
    //     ampo += -0.5*h, "S+" , j ;
    //     ampo += -0.5*h, "S-" , j ;
    //     ampo += -0.5*h, "S+" , j+1 ;
    //     ampo += -0.5*h, "S-" , j+1 ;
    // }

    // // chain 2 PBC
    //     ampo += J, "Sz", N, "Sz", 2 ;
    
    // // Hopping 
    // for (int j = 1; j <= N-1 ; j += 3) 
    // {
    //    ampo += g * 0.5, "S+", j, "Sz", j+1 ; 
    //    ampo += g * 0.5, "S-", j, "Sz", j+1 ; 
    // }
    // for (int j = 1; j <= N-4 ; j += 3) 
    // {
    //   ampo += -g, "Sz", j, "Sz", j+3 ,"Sz" , j+2 ;
    // }
    
    // // Hopping PBC
    //   ampo += -g, "Sz", N-2, "Sz", 1 ,"Sz" , N ;

    H = toMPO(ampo);

    // Several tries to the occupation of initial state
    //
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%6 == 1)             state.set(i,"Up");
        else if(i%6 == 4)        state.set(i,"Up");
        else if(i%6 == 2)        state.set(i,"Up");
        else if(i%6 == 5)        state.set(i,"Dn");
        else if(i%6 == 3)        state.set(i,"Up");
        else if(i%6 == 0)        state.set(i,"Dn");
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
    
//Given an MPS called "psi",
//nd some particular bond "b" (1 <= b < length(psi))
//across which we want to compute the von Neumann entanglement

for(int b=3 ; b<=N-1 ; b +=3 )
{
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
printfln("SvN = ", SvN);
}
    return 0;
}