#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    int N = 4;
    MPO H;
    MPS psi0;
    MPO H1 ;
    MPO n_j ;
    //
    // Initialize the site degrees of freedom
    // Setting "ConserveQNs=",true makes the indices
    // carry Sz quantum numbers and will lead to 
    // block-sparse MPO and MPS tensors
    //

    // kitaev model 在 lieb lattice 上面 ，由于 kitaev model p-wave 超导，自旋同向。
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    auto sites = Fermion(N,{"ConserveQNs=", false}); 


    auto t=2;
    auto D=2;
    auto V=2;
    auto J=2;

    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-2 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+2, "C", j;
            ampo += -t, "Cdag", j,   "C", j+2;
        
    }
    // for (int j = 1; j <= N-1 ; j += 1) //electron pairing
    //     {
    //        ampo += D, "Adag", j+1, "Adag", j;
    //        ampo += D, "A", j,   "A", j+1;
    // }
    for (int j = 1; j <= N-2 ; j += 1) //electron pairing
        {
           ampo += D, "Cdag", j+2, "Cdag", j;
           ampo += D, "C", j,   "C", j+2;
    }

    for (int j = 1; j <= N ; j += 1) //onsite V 
    {
        ampo += 0.5*V, "Cdag", j, "C", j;
        ampo += -0.5*V, "C", j, "Cdag", j;
    }

    // for (int j = 2; j <= N ; j += 2) //onsite V
    // {
    //     ampo += -0.5*V, "Cdag", j, "C", j;
    //     ampo +=  0.5*V, "C", j, "Cdag", j;
    // }

    ampo += -J, "Cdag", 3, "C", 2;
    ampo += -J, "Cdag", 2, "C", 3;

    auto state = InitState(sites);



    for(auto i : range1(N))
        {
        if(i%2 == 1) state.set(i,"Emp");
        else         state.set(i,"Occ");
        }
    psi0 = MPS(state);
    H = toMPO(ampo);

    auto ampo1 = AutoMPO(sites);
    for (int j = 1; j <= N-2 ; j += 1) //electron hopping 
        {
            ampo1 += -t, "Cdag", j+2, "C", j;
            ampo1 += -t, "Cdag", j,   "C", j+2;
        
    }
    // for (int j = 1; j <= N-1 ; j += 1) //electron pairing
    //     {
    //        ampo += D, "Adag", j+1, "Adag", j;
    //        ampo += D, "A", j,   "A", j+1;
    // }
    for (int j = 1; j <= N-2 ; j += 1) //electron pairing
        {
           ampo1 += D, "Cdag", j+2, "Cdag", j;
           ampo1 += D, "C", j,   "C", j+2;
    }

    for (int j = 1; j <= N ; j += 1) //onsite V 
    {
        ampo1 += 0.5*V, "Cdag", j, "C", j;
        ampo1 += -0.5*V, "C", j, "Cdag", j;
    }

    // for (int j = 2; j <= N ; j += 2) //onsite V
    // {
    //     ampo += -0.5*V, "Cdag", j, "C", j;
    //     ampo +=  0.5*V, "C", j, "Cdag", j;
    // }

    ampo1 += -J, "Cdag", 3, "C", 2;


    H1 = toMPO(ampo1);

    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10,20,100,100,200 ;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
    auto [energy1,psi1] = dmrg(H1,psi0,sweeps,"Quiet");
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\npsi = %.10f",psi );
    printfln("\npsi1 = %.10f",psi1 );
    // psi.position(3) ;
    // printfln("\nUsing inner = %.10f",psi );
    // auto Sz_2 = op(sites,"C",2);
    // printfln("\nUsing inner = %.10f",Sz_2 );
    auto b = inner(psi,psi1) ; 
    printfln("\nb= %.10f",b);
    


    // itensor 输出的psi 是一个MPS 但是可以通过计算 n_i 的期望值得到近似特征向量。
    for(int j = 1; j <= N; ++j)
        {
        //re-gauge psi to get ready to measure at position j
        auto ampo2 = AutoMPO(sites) ;
        ampo2 += 1, "Cdag", j, "C", j;
        n_j = toMPO(ampo2);
        auto wave = inner(psi , n_j , psi) ;
        printfln("Site %d occupation: %.5f", j, wave);
    }
    
    return 0;
}