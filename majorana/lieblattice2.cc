#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    // N must be 3Z
    int N = 15;
    MPO H;
    MPS psi0;
    MPO Hpbc ;
    MPO n_j ;
    MPO G1;
    MPO G2;
    MPO G3;
    MPO G4;
    //
    // Initialize the site degrees of freedom
    // Setting "ConserveQNs=",true makes the indices
    // carry Sz quantum numbers and will lead to 
    // block-sparse MPO and MPS tensors
    //

    // kitaev model 在 lieb lattice 上面 ，由于 kitaev model p-wave 超导，自旋同向。
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    auto sites = Fermion(N,{"ConserveQNs=", false}); 


    auto t= 0.03;
    auto J= 1;
    auto V= 0.2;
    auto mu = 0; 
    auto D =0.6*t;
    auto U = 0.01 ;
    // auto t=1;
    // auto J=1;
    // auto V=1;
    // auto mu = 0; 
    // auto D =1 ;

    // Open boundary condition
    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-3 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+3, "C", j;
            ampo += -t, "Cdag", j,   "C", j+3;
    }
    for(int j = 1; j <= N-5; j += 3) // electron hopping over orbitals
        {
            ampo += -J, "Cdag", j+1, "C", j;
            ampo += -J, "Cdag", j,   "C", j+1;
            ampo += -J, "Cdag", j+2, "C", j;
            ampo += -J, "Cdag", j,   "C", j+2;
            ampo += -J, "Cdag", j+4, "C", j;
            ampo += -J, "Cdag", j,   "C", j+4;
            ampo += -J, "Cdag", j+5, "C", j;
            ampo += -J, "Cdag", j,   "C", j+5;
    }
        
    for (int j = 2; j <= N-1 ; j += 3) //onsite V 这里几项我尝试配做particle-hole对称的形式，去掉了一个常数项
        {
            ampo += V, "Cdag", j, "C", j;
          //  ampo += -0.5*V, "C", j, "Cdag", j;
    }

    for (int j = 3; j <= N ; j += 3) //onsite V
        {
            ampo += -V, "Cdag", j, "C", j;
          //  ampo +=  0.5*V, "C", j, "Cdag", j;
    }

    for (int j = 1; j <= N ; j += 1) //chemical potential
        {
            // ampo += -0.5*mu-0.0001, "Cdag", j, "C", j;
            // ampo += 0.5*mu-0.0001, "C", j, "Cdag", j;
            ampo += -0.5*mu, "Cdag", j, "C", j;
            ampo += 0.5*mu, "C", j, "Cdag", j;
    }

    for (int j = 1; j <= N-3 ; j += 1) //electron pairing
        {
            ampo += D, "Cdag", j+3, "Cdag", j;
            ampo += D, "C", j,   "C", j+3;
    }


    // interaction
    for (int j = 1; j <= N-3 ; j += 1) 
        {
            ampo += U , "Cdag", j, "C", j , "Cdag" , j+3 , "C" , j+3 ;
    }
    
    ampo += -J, "Cdag", N-1, "C", N-2;
    ampo += -J, "Cdag", N-2,   "C", N-1;
    ampo += -J, "Cdag", N, "C", N-2;
    ampo += -J, "Cdag", N-2,   "C", N;

    H = toMPO(ampo);
    // Several tries to the occupation of initial state
    //
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%3 == 0) state.set(i,"Occ");
        else         state.set(i,"Emp");
        }
    psi0 = MPS(state);
    
    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>

    // sweep
    auto sweeps = Sweeps(50);
    sweeps.maxdim() = 10,10,10,10,20,20,20,20,40,40,40,40,100,100,100,100,200,200,200,200,400,400,400,400,600,600,600,600;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = -1E-3,-1E-4,-1E-5,-1E-6,-1E-7,-1E-9,0.0;
    //sweeps.noise() = 1E-1,1E-1,1E-2,1E-2,1E-2,1E-2,1E-2,1E-2,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-4,1E-4,1E-4,1E-4,1E-4,1E-4,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-9,1E-9,1E-9,1E-9,1E-9,1E-9,1E-10,1E-10,1E-10,1E-10,1E-10,1E-10,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,0.0;
    //println(sweeps);


    // Begin the DMRG calculation
   auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
   printfln("Initial energy = %.5f", inner(psi0,H,psi0) );
   printfln("\nGround State Energy = %.10f",energy);
   

    // calculating the pair-correlation function based on tensor contraction 
    for(int j = 2; j < N-1; ++j)
        {
        auto C1 = op(sites,"Cdag",1);
        auto Cj = op(sites,"Cdag",j);
        psi.position(1); 
        auto psidag = dag(psi);
        psidag.prime("Link");
        //auto li_1 = leftLinkIndex(psi,1);
        //auto C = prime(psi(1),li_1)*C1;
        auto C = psi(1)*C1;
        C *= prime(psidag(1),"Site");
        for(int k = 2; k < j; ++k)
            {
            C *= psi(k);
            C *= psidag(k);
            }
        auto lj = rightLinkIndex(psi,j);
        C *= prime(psi(j),lj)*Cj;
        C *= prime(psidag(j),"Site");
        auto result = elt(C); 
        print(result) ;
        }

    
    //printfln("\nGround State Energy = %.10f",result);



   return 0 ;
    }