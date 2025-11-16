#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    int N = 180;
    MPO H;
    MPS psi0;
    MPO H1 ;
    MPO n_j ;
    int mid =N/4 ;
    //
    // Initialize the site degrees of freedom
    // Setting "ConserveQNs=",true makes the indices
    // carry Sz quantum numbers and will lead to 
    // block-sparse MPO and MPS tensors
    //

    // kitaev model 在 lieb lattice 上面 ，由于 kitaev model p-wave 超导，自旋同向。
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    auto sites = Fermion(N,{"ConserveQNs=", false}); 


    auto t=0.06;
    auto D=0.03;

    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-3 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+3, "C", j;
            ampo += -t, "Cdag", j,   "C", j+3;
        
    }
    // for (int j = 1; j <= N-1 ; j += 1) //electron pairing
    //     {
    //        ampo += D, "Adag", j+1, "Adag", j;
    //        ampo += D, "A", j,   "A", j+1;
    // }
    for (int j = 1; j <= N-3 ; j += 1) //electron pairing
        {
           ampo += D, "Cdag", j+3, "Cdag", j;
           ampo += D, "C", j,   "C", j+3;
    }

    // for (int j = 1; j <= N ; j += 1) //onsite V 
    // {
    //     ampo += 0.5*V, "Cdag", j, "C", j;
    //     ampo += -0.5*V, "C", j, "Cdag", j;
    // }

    // for (int j = 2; j <= N ; j += 2) //onsite V
    // {
    //     ampo += -0.5*V, "Cdag", j, "C", j;
    //     ampo +=  0.5*V, "C", j, "Cdag", j;
    // }

    // ampo += -J, "Cdag", 3, "C", 2;
    // ampo += -J, "Cdag", 2, "C", 3;

    auto state = InitState(sites);



    for(auto i : range1(N))
        {
        if(i%2 == 1) state.set(i,"Occ");
        else         state.set(i,"Occ");
        }
    psi0 = MPS(state);
    H = toMPO(ampo);

    // auto ampo1 = AutoMPO(sites);
    // for (int j = 1; j <= N-1 ; j += 1) //electron hopping 
    //     {
    //         ampo1 += -t, "Cdag", j+1, "C", j;
    //         ampo1 += -t, "Cdag", j,   "C", j+1;
        
    // }
    // // for (int j = 1; j <= N-1 ; j += 1) //electron pairing
    // //     {
    // //        ampo += D, "Adag", j+1, "Adag", j;
    // //        ampo += D, "A", j,   "A", j+1;
    // // }
    // for (int j = 1; j <= N-1 ; j += 1) //electron pairing
    //     {
    //        ampo1 += D, "Cdag", j+1, "Cdag", j;
    //        ampo1 += D, "C", j,   "C", j+1;
    // }

    // // for (int j = 1; j <= N ; j += 1) //onsite V 
    // // {
    // //     ampo1 += 0.5*V, "Cdag", j, "C", j;
    // //     ampo1 += -0.5*V, "C", j, "Cdag", j;
    // // }

    // // for (int j = 2; j <= N ; j += 2) //onsite V
    // // {
    // //     ampo += -0.5*V, "Cdag", j, "C", j;
    // //     ampo +=  0.5*V, "C", j, "Cdag", j;
    // // }

    // // ampo1 += -J, "Cdag", 3, "C", 2;
    // ampo1 += -t, "Cdag", N, "C", 1;
    // ampo1 += -t, "Cdag", 1, "C", N;
    // ampo1 += D, "Cdag", N, "Cdag", 1;
    // ampo1 += D, "C", 1, "C", N;

    // H1 = toMPO(ampo1);

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
    auto sweeps = Sweeps(30);
    sweeps.maxdim() = 10,20,100,100,200,400 ;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-5,1E-6,1E-7,1E-8,0.0;
    //println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
    //auto [energy1,psi1] = dmrg(H1,psi0,sweeps,"Quiet");
    printfln("\nGround State Energy = %.10f",energy);
    //printfln("\nGround State Energy = %.10f",energy1);
    // printfln("\npsi = %.10f",psi );
    // printfln("\npsi1 = %.10f",psi1 );
    // psi.position(3) ;
    // printfln("\nUsing inner = %.10f",psi );
    // auto Sz_2 = op(sites,"C",2);
    // printfln("\nUsing inner = %.10f",Sz_2 );
    
    std::ofstream outfile5("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try8");
    for(int b = 3 ; b < N-1; b +=1 )
        {
        psi.position(b); 

        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U,S,V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);

        Real SvN = 0.;
        for(auto n : range1(dim(u)))
            {
            auto Sn = elt(S,n,n);
            auto p = sqr(Sn);
            if(p > 1E-12) SvN += -p*log(p);
            }
        printfln("Across bond b=%d, SvN = %.10f",b,SvN);
        outfile5 << SvN << std::endl;
        }
    outfile5.close();

    // // itensor 输出的psi 是一个MPS 但是可以通过计算 n_i 的期望值得到近似特征向量。
    // for(int j = 1; j <= N; ++j)
    //     {
    //     //re-gauge psi to get ready to measure at position j
    //     auto ampo2 = AutoMPO(sites) ;
    //     ampo2 += 1, "Cdag", j, "C", j;
    //     n_j = toMPO(ampo2);
    //     auto wave = inner(psi , n_j , psi) ;
    //     printfln("Site %d occupation: %.5f", j, wave);
    // }
    
    // for(int j = 1; j <= N; ++j)
    //     {
    //     //re-gauge psi to get ready to measure at position j
    //     auto ampo2 = AutoMPO(sites) ;
    //     ampo2 += 1, "Cdag", j, "C", j;
    //     n_j = toMPO(ampo2);
    //     auto wave1 = inner(psi , n_j , psi) ;
    //     auto wave2 = inner(psi1 , n_j , psi) ;
    //     auto wave3 = inner(psi , n_j , psi1) ;
    //     auto wave4 = inner(psi1 , n_j , psi1) ;
    //     auto b     = inner(psi1,psi) ;
    //     auto majorana = wave1 - b*wave2 - b*wave3 + b*b*wave4 ;
    //     printfln("Site %d occupation: %.5f", j, majorana);
    // }

    std::ofstream outfile9("/mnt/d/OneDrive/programs/vscodepython/marjorana/output7");
    for(int j = mid+3; j < N-1; j+=3)
        {
        auto Adag3 = op(sites,"Adag",mid);
        auto A3 = op(sites,"A",mid);
        auto Aj = op(sites,"A",j);
        auto Adagj = op(sites,"Adag",j);


        // guage psi is a must for contracting left side
        psi.position(mid) ;
        auto psidag = dag(psi);
        psidag.prime();
        auto li_1 = leftLinkIndex(psi,mid);

        //constructing majorana operator using spinless fermion basis
        auto Adag3Aj = prime(psi(mid),li_1)*Adag3*psidag(mid);
        auto  A3Aj = prime(psi(mid),li_1)*A3*psidag(mid);
        auto Adag3Adagj = prime(psi(mid),li_1)*Adag3*psidag(mid);
        auto  A3Adagj = prime(psi(mid),li_1)*A3*psidag(mid);

        for(int k = mid+1; k < j; ++k)
            {
            Adag3Aj *= psi(k);
            Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
            Adag3Aj *= psidag(k);
            A3Aj *= psi(k);
            A3Aj *= op(sites,"F",k); //Jordan-Wigner string
            A3Aj *= psidag(k);
            Adag3Adagj *= psi(k);
            Adag3Adagj *= op(sites,"F",k); //Jordan-Wigner string
            Adag3Adagj *= psidag(k);
            A3Adagj *= psi(k);
            A3Adagj *= op(sites,"F",k); //Jordan-Wigner string
            A3Adagj *= psidag(k);
            }
        auto lj = rightLinkIndex(psi,j);

        Adag3Aj  *= prime(psi(j),lj);
        Adag3Aj  *= Aj;
        Adag3Aj  *= psidag(j);
        A3Aj  *= prime(psi(j),lj);
        A3Aj  *= Aj;
        A3Aj  *= psidag(j);
        Adag3Adagj  *= prime(psi(j),lj);
        Adag3Adagj  *= Adagj;
        Adag3Adagj  *= psidag(j);
        A3Adagj  *= prime(psi(j),lj);
        A3Adagj  *= Adagj;
        A3Adagj  *= psidag(j);

        //consider JW-transformaton 
        //origin : c3cj - c3cdagj + cdag3cj - cdag3cdagj
        //now:    -a3aj + a3adagj + adag3aj - adag3agdaj 
        //auto result =   elt(A3Adagj) ;//- elt(A3Aj) - elt(Adag3Adagj) + elt(Adag3Aj);
        auto result = - elt(A3Aj) ;
        outfile9  << result << std::endl;
        }
    outfile9.close();

    std::ofstream outfile10("/mnt/d/OneDrive/programs/vscodepython/marjorana/output10");
    for(int j = mid+3; j < N-1; j+=3)
        {
        auto Adag3 = op(sites,"Adag",mid);
        auto A3 = op(sites,"A",mid);
        auto Aj = op(sites,"A",j);
        auto Adagj = op(sites,"Adag",j);


        // guage psi is a must for contracting left side
        psi.position(mid) ;
        auto psidag = dag(psi);
        psidag.prime();
        auto li_1 = leftLinkIndex(psi,mid);

        //constructing majorana operator using spinless fermion basis
        auto Adag3Aj = prime(psi(mid),li_1)*Adag3*psidag(mid);
        auto  A3Aj = prime(psi(mid),li_1)*A3*psidag(mid);
        auto Adag3Adagj = prime(psi(mid),li_1)*Adag3*psidag(mid);
        auto  A3Adagj = prime(psi(mid),li_1)*A3*psidag(mid);

        for(int k = mid+1; k < j; ++k)
            {
            Adag3Aj *= psi(k);
            Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
            Adag3Aj *= psidag(k);
            A3Aj *= psi(k);
            A3Aj *= op(sites,"F",k); //Jordan-Wigner string
            A3Aj *= psidag(k);
            Adag3Adagj *= psi(k);
            Adag3Adagj *= op(sites,"F",k); //Jordan-Wigner string
            Adag3Adagj *= psidag(k);
            A3Adagj *= psi(k);
            A3Adagj *= op(sites,"F",k); //Jordan-Wigner string
            A3Adagj *= psidag(k);
            }
        auto lj = rightLinkIndex(psi,j);

        Adag3Aj  *= prime(psi(j),lj);
        Adag3Aj  *= Aj;
        Adag3Aj  *= psidag(j);
        A3Aj  *= prime(psi(j),lj);
        A3Aj  *= Aj;
        A3Aj  *= psidag(j);
        Adag3Adagj  *= prime(psi(j),lj);
        Adag3Adagj  *= Adagj;
        Adag3Adagj  *= psidag(j);
        A3Adagj  *= prime(psi(j),lj);
        A3Adagj  *= Adagj;
        A3Adagj  *= psidag(j);

        //consider JW-transformaton 
        //origin : c3cj - c3cdagj + cdag3cj - cdag3cdagj
        //now:    -a3aj + a3adagj + adag3aj - adag3agdaj 
        auto result =   elt(A3Adagj) ;//- elt(A3Aj) - elt(Adag3Adagj) + elt(Adag3Aj);
        //auto result = - elt(A3Aj) ;
        outfile10  << result << std::endl;
        }
    outfile10.close();

    return 0;
}