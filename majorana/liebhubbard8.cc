#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    // N must be 3Z
    int N = 450 ;
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

    // attractive hubbard model
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    auto sites = Fermion(N,{"ConserveQNs=", false}); 

    auto t= 0.03;
    auto J= 1;
    auto V= 0.1;
    auto mu = 0; 
    auto D = 0.6*t ;

    // Open boundary condition
    auto ampo = AutoMPO(sites);

    // the left half
    for (int j = 1; j <= N/2 -3 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+3, "C", j;
            ampo += -t, "Cdag", j,   "C", j+3;
        }
    for(int j = 3; j <= N/2-3; j += 3) // electron hopping over orbitals
        {
            ampo += -J, "Cdag", j-2, "C", j;
            ampo += -J, "Cdag", j,   "C", j-2;
            ampo += -J, "Cdag", j+1, "C", j;
            ampo += -J, "Cdag", j,   "C", j+1;
            ampo += -J, "Cdag", j-1, "C", j;
            ampo += -J, "Cdag", j,   "C", j-1;
            ampo += -J, "Cdag", j+2, "C", j;
            ampo += -J, "Cdag", j,   "C", j+2;
    }
        
    for (int j = 1; j <= N/2-2 ; j += 3) //onsite V 
        {
            ampo += V, "Cdag", j, "C", j;
    }

    for (int j = 2; j <= N/2 -1 ; j += 3) //onsite V
        {
            ampo += -V, "Cdag", j, "C", j;
    }

    // electron pairing
    for (int j = 1; j <= N/2-3 ; j += 1) //electron pairing
        {
            ampo += D, "Cdag", j+3, "Cdag", j;
            ampo += D, "C", j,   "C", j+3;
    }
   
    // right half
    for (int j = N/2 + 1 ; j <= N -3 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+3, "C", j;
            ampo += -t, "Cdag", j,   "C", j+3;
        }
    for(int j = N/2 + 1; j <= N-5; j += 3) // electron hopping over orbitals
        {
            ampo += -J, "Cdag", j+2, "C", j;
            ampo += -J, "Cdag", j,   "C", j+2;
            ampo += -J, "Cdag", j+1, "C", j;
            ampo += -J, "Cdag", j,   "C", j+1;
            ampo += -J, "Cdag", j+4, "C", j;
            ampo += -J, "Cdag", j,   "C", j+4;
            ampo += -J, "Cdag", j+5, "C", j;
            ampo += -J, "Cdag", j,   "C", j+5;
    }

       
    for (int j = N/2 +3; j <= N ; j += 3) //onsite V 
        {
            ampo += V, "Cdag", j, "C", j;
    }

    for (int j = N/2 +2; j <= N-1 ; j += 3) //onsite V
        {
            ampo += -V, "Cdag", j, "C", j;
    }

    // electron pairing
    for (int j = N/2 + 1 ; j <= N -3 ; j += 1) //electron pairing
        {
            ampo += D, "Cdag", j+3, "Cdag", j;
            ampo += D, "C", j,   "C", j+3;
    }

    // middle
    ampo += -t, "Cdag", N/2, "C", N/2+1;
    ampo += -t, "Cdag", N/2+1,   "C", N/2;
    ampo += -t, "Cdag", N/2-1, "C", N/2+2;
    ampo += -t, "Cdag", N/2+2,   "C", N/2-1;
    ampo += -t, "Cdag", N/2-2, "C", N/2+3;
    ampo += -t, "Cdag",N/2+3,   "C", N/2-2;
    ampo += D, "Cdag", N/2 +1 , "Cdag", N/2;
    ampo += D, "C", N/2,   "C", N/2 +1;
    ampo += D, "Cdag", N/2+2, "Cdag", N/2-1;
    ampo += D, "C", N/2-1,   "C", N/2+2;
    ampo += D, "Cdag", N/2+3, "Cdag", N/2-2;
    ampo += D, "C", N/2-2,   "C", N/2+3;
    ampo += -J, "Cdag", N/2, "C", N/2-2;
    ampo += -J, "Cdag", N/2-2,   "C", N/2;
    ampo += -J, "Cdag", N/2, "C", N/2-1;
    ampo += -J, "Cdag", N/2-1,   "C", N/2;
    ampo += -J, "Cdag", N/2, "C", N/2+3;
    ampo += -J, "Cdag", N/2+3,   "C", N/2;
    ampo += -J, "Cdag", N/2, "C", N/2+2;
    ampo += -J, "Cdag", N/2+2,   "C", N/2;

    //LAST
    ampo += -J, "Cdag", N-2, "C", N;
    ampo += -J, "Cdag", N,   "C", N-2;
    ampo += -J, "Cdag", N-2, "C", N-1;
    ampo += -J, "Cdag", N-1,   "C", N-2;


    H = toMPO(ampo);
    // Several tries to the occupation of initial state
    //
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%3 == 3) state.set(i,"Occ");
        else         state.set(i,"Emp");
        }
    psi0 = MPS(state);
    
    // sweep
    auto sweeps = Sweeps(30);
    sweeps.maxdim() = 10,20,40,100,200,400,800,1600,3200,6400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,1E-9,1E-10,0.0;

    // sweep
    // auto sweeps = Sweeps(150);
    // sweeps.maxdim() = 10,10,10,10,20,20,20,20,40,40,40,40,100,100,100,100,200,200,200,200,400,400,400,400,600,600,600,600;
    // sweeps.cutoff() = 1E-10;
    // sweeps.niter() = 2;
    // sweeps.noise() = -1E-3,-1E-4,-1E-5,-1E-6,-1E-7,-1E-9,0.0;
    // //sweeps.noise() = 1E-1,1E-1,1E-2,1E-2,1E-2,1E-2,1E-2,1E-2,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-4,1E-4,1E-4,1E-4,1E-4,1E-4,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-9,1E-9,1E-9,1E-9,1E-9,1E-9,1E-10,1E-10,1E-10,1E-10,1E-10,1E-10,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,0.0;
    //println(sweeps);


    // Begin the DMRG calculation
   auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
   printfln("Initial energy = %.5f", inner(psi0,H,psi0) );
   printfln("\nGround State Energy = %.10f",energy);

    // store
    std::string filename = "outputpsi/asymmetry_order2_sweep10_N_" 
                        + std::to_string(N) + "_t_" + std::to_string(t) 
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V) 
                        + "_D_" + std::to_string(D) ;
    std::string filename2 = "outputpsi/sites_order2_sweep10_N_" 
                        + std::to_string(N) + "_t_" + std::to_string(t) 
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V) 
                        + "_D_" + std::to_string(D) ;
                        
    writeToFile(filename2,sites);
    writeToFile(filename, psi);
    //calulating on-site ni
    // for(int j = 1; j < N; ++j) 
    //     {
    //     psi.position(j);
    //     auto ket = psi(j);
    //     auto bra = dag(prime(ket,"Site"));
    //     auto nj = op(sites,"Cdag",j)*op(sites,"C",j);
    //     auto occupation = elt(bra*nj*ket);
    //     print(occupation);
    //     }
    

    // measuring majorana with poor efficiency , used only for testing
    // std::ofstream outfile("output");
    // for(int j = 1; j <= N; ++j)
    //     {
    //     // //re-gauge psi to get ready to measure at position j
    //     // auto ampo3 = AutoMPO(sites) ;
    //     // ampo3 += 1, "C", 3, "C", j;
    //     // G1 = toMPO(ampo3);
    //     // auto ampo4 = AutoMPO(sites) ;
    //     // ampo4 += 1, "Cdag", 3, "C", j;
    //     // G2 = toMPO(ampo4);
    //     auto ampo5 = AutoMPO(sites) ;
    //     ampo5 += -1, "C", j, "C", N;
    //     G3 = toMPO(ampo5);
    //     auto ampo6 = AutoMPO(sites) ;
    //     ampo6 += 1, "C", j, "Cdag", N;
    //     G4 = toMPO(ampo6);
    //     // auto wave1 = inner(psi , G1 , psi) ;
    //     // auto wave2 = inner(psi , G2 , psi) ;
    //     auto wave3 = inner(psi , G3 , psi) ;
    //     auto wave4 = inner(psi , G4 , psi) ;
    //     auto majorana= wave3 + wave4 ;//+ wave3 + wave4 ;
    //     printfln("Site %d occupation: %.5f", j, majorana);
    //     outfile << majorana << std::endl;
    //     }
    //     //print(result) ; 
    // outfile.close();
    

    //calcualting left majorana <c_j gamma_3>
    std::ofstream outfile3("/mnt/d/OneDrive/programs/vscodepython/marjorana/output3");
    for(int j = 6; j < (N+1)/2; j +=3)
        {
        auto Adag3 = op(sites,"Adag",3);
        auto A3 = op(sites,"A",3);
        auto Aj = op(sites,"A",j);
        
        // guage psi is a must for contracting left side
        psi.position(3) ;
        auto psidag = dag(psi);
        psidag.prime();
        auto li_1 = leftLinkIndex(psi,3);
        auto Adag3Aj = prime(psi(3),li_1)*Adag3*psidag(3);
        auto  A3Aj = prime(psi(3),li_1)*A3*psidag(3);
        for(int k = 4; k < j; ++k)
            {
            Adag3Aj *= psi(k);
            Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
            Adag3Aj *= psidag(k);
            A3Aj *= psi(k);
            A3Aj *= op(sites,"F",k); //Jordan-Wigner string
            A3Aj *= psidag(k);
            }
        auto lj = rightLinkIndex(psi,j);
        Adag3Aj  *= prime(psi(j),lj);
        Adag3Aj  *= Aj;
        Adag3Aj  *= psidag(j);
        A3Aj  *= prime(psi(j),lj);
        A3Aj  *= Aj;
        A3Aj  *= psidag(j);

        //consider JW-transformaton the second sign is minus
        auto result = elt(Adag3Aj) - elt(A3Aj);
        outfile3 << result << std::endl;
        }

    for(int j =  (N+1)/2 +1; j < N-2; j +=3)
        {
        auto Adag3 = op(sites,"Adag",3);
        auto A3 = op(sites,"A",3);
        auto Aj = op(sites,"A",j);
        
        // guage psi is a must for contracting left side
        psi.position(3) ;
        auto psidag = dag(psi);
        psidag.prime();
        auto li_1 = leftLinkIndex(psi,3);
        auto Adag3Aj = prime(psi(3),li_1)*Adag3*psidag(3);
        auto  A3Aj = prime(psi(3),li_1)*A3*psidag(3);
        for(int k = 4; k < j; ++k)
            {
            Adag3Aj *= psi(k);
            Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
            Adag3Aj *= psidag(k);
            A3Aj *= psi(k);
            A3Aj *= op(sites,"F",k); //Jordan-Wigner string
            A3Aj *= psidag(k);
            }
        auto lj = rightLinkIndex(psi,j);
        Adag3Aj  *= prime(psi(j),lj);
        Adag3Aj  *= Aj;
        Adag3Aj  *= psidag(j);
        A3Aj  *= prime(psi(j),lj);
        A3Aj  *= Aj;
        A3Aj  *= psidag(j);

        //consider JW-transformaton the second sign is minus
        auto result = elt(Adag3Aj) - elt(A3Aj);
        outfile3 << result << std::endl;
        }
    outfile3.close();
     
    // //calcualting right majorana <c_j*i*gamma_N>
    // std::ofstream outfile4("/mnt/d/OneDrive/programs/vscodepython/marjorana/output4");
    // for(int j = 4; j < N-1; ++j)
    //     {
    //     auto AdagN = op(sites,"Adag",N);
    //     auto AN = op(sites,"A",N);
    //     auto Aj = op(sites,"A",j);

    //     psi.position(j) ;
    //     auto psidag = dag(psi);
    //     psidag.prime();
    //     auto li_j = leftLinkIndex(psi,j);
    //     auto AjAdagN = prime(psi(j),li_j)*Aj*psidag(j);
    //     auto  AjAN = prime(psi(j),li_j)*Aj*psidag(j);
    //     for(int k = j+1 ; k < N; ++k)
    //         {
    //         AjAdagN *= psi(k);
    //         AjAdagN *= op(sites,"F",k); //Jordan-Wigner string
    //         AjAdagN *= psidag(k);
    //         AjAN *= psi(k);
    //         AjAN *= op(sites,"F",k); //Jordan-Wigner string
    //         AjAN *= psidag(k);
    //         }
    //     AjAdagN  *= psi(N) ;
    //     AjAdagN  *= AdagN;
    //     AjAdagN  *= psidag(N);
    //     AjAN  *= psi(N) ;
    //     AjAN  *= AN;
    //     AjAN  *= psidag(N);
       
    //     //consider JW-transformaton the second sign is plus ,first sign is minus
    //     auto result = - elt(AjAdagN) + elt(AjAN);
    //     outfile4  << "\t" << result << std::endl;
    //     }
    // outfile4.close();

    //calculate entanglment entropy
    std::ofstream outfile5("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try1");
    for(int b = 3 ; b <= N/2 -3 ; b +=3 )
    //for(int b = 3 ; b <= N-2 ; b +=3 )
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
    
    for(int b = N/2 + 3  ; b <= N-3 ; b +=3 )
    //for(int b = N/2 + 4  ; b <= N-3 ; b +=3 )
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

    std::ofstream outfile6("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try2");
    for(int b = 4 ; b <= N/2 - 2 ; b +=3 )
    //for(int b = 2 ; b <= N/2 - 2 ; b +=3 )
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
        outfile6 << SvN << std::endl;
        }
    
    for(int b = N/2 + 2 ; b <= N-4 ; b +=3 )
    //for(int b = N/2 + 3 ; b <= N-4 ; b +=3 )
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
        outfile6 << SvN << std::endl;
        }
    
    outfile6.close();

    std::ofstream outfile7("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try3");
    for(int b = 2 ; b <= N/2 - 1 ; b +=3 )
    //for(int b = 3 ; b <= N/2 - 1 ; b +=3 )
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
        outfile7 << SvN << std::endl;
        }
    
    for(int b = N/2 +1 ; b <= N-2 ; b +=3 )
    //for(int b = N/2 +2 ; b <= N-2 ; b +=3 )
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
        outfile7 << SvN << std::endl;
        }
    
    outfile7.close();





   return 0 ;
    }