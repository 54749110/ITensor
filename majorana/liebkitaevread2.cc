#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    
    int N = 150 ;
    auto sites = Fermion(N,{"ConserveQNs=", false});
    readFromFile("outputpsi/sites_N_try2_150_t_0.030000_J_1_V_0.100000_D_0.018000",sites);
    MPS psi(sites);
    readFromFile<MPS>("outputpsi/N_try2_150_t_0.030000_J_1_V_0.100000_D_0.018000", psi);



    // std::ofstream outfile3("/mnt/d/OneDrive/programs/vscodepython/marjorana/output3");
    // for(int j = 6; j < N-1; j +=3)
    //     {
    //     auto Adag3 = op(sites,"Adag",3);
    //     auto A3 = op(sites,"A",3);
    //     auto Aj = op(sites,"A",j);
        
    //     // guage psi is a must for contracting left side
    //     psi.position(3) ;
    //     auto psidag = dag(psi);
    //     psidag.prime();
    //     auto li_1 = leftLinkIndex(psi,3);
    //     auto Adag3Aj = prime(psi(3),li_1)*Adag3*psidag(3);
    //     auto  A3Aj = prime(psi(3),li_1)*A3*psidag(3);
    //     for(int k = 4; k < j; ++k)
    //         {
    //         Adag3Aj *= psi(k);
    //         Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
    //         Adag3Aj *= psidag(k);
    //         A3Aj *= psi(k);
    //         A3Aj *= op(sites,"F",k); //Jordan-Wigner string
    //         A3Aj *= psidag(k);
    //         }
    //     auto lj = rightLinkIndex(psi,j);
    //     Adag3Aj  *= prime(psi(j),lj);
    //     Adag3Aj  *= Aj;
    //     Adag3Aj  *= psidag(j);
    //     A3Aj  *= prime(psi(j),lj);
    //     A3Aj  *= Aj;
    //     A3Aj  *= psidag(j);

    //     //consider JW-transformaton the second sign is minus
    //     auto result = elt(Adag3Aj) - elt(A3Aj);
    //     outfile3  << result << std::endl;
    //     }
    // outfile3.close();

 //calculate entanglment entropy
    std::ofstream outfile5("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try1");
    //for(int b = 3 ; b <= N/2 -3 ; b +=3 )
    for(int b = 3 ; b <= N-2 ; b +=3 )
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
    
    // for(int b = N/2 + 3  ; b <= N-3 ; b +=3 )
    // //for(int b = N/2 + 4  ; b <= N-3 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile5 << SvN << std::endl;
    //     }
    
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


    // even symmetry
    //calculate entanglment entropy
    // std::ofstream outfile5("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try4");
    // for(int b = 3 ; b <= (N+1)/2 - 2 ; b +=3 )
    // // for(int b = 4 ; b <= (N+1)/2 - 2 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile5 << SvN << std::endl;
    //     }
    
    // for(int b = (N+1)/2 + 2 ; b <= N-2 ; b +=3 )
    // // for(int b = (N+1)/2 + 3 ; b <= N-2 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile5 << SvN << std::endl;
    //     }
    
    // outfile5.close();

    // std::ofstream outfile6("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try5");
    // for(int b = 4 ; b <= (N+1)/2 - 2 ; b +=3 )
    // // for(int b = 2 ; b <= (N+1)/2 - 2 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile6 << SvN << std::endl;
    //     }
    
    // for(int b = (N+1)/2 + 1 ; b <= N-2 ; b +=3 )
    // // for(int b = (N+1)/2 + 2 ; b <= N-2 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile6 << SvN << std::endl;
    //     }
    
    // outfile6.close();

    // std::ofstream outfile7("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try6");
    // for(int b = 2 ; b <= (N+1)/2 - 2 ; b +=3 )
    // // for(int b = 3 ; b <= (N+1)/2 - 2 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile7 << SvN << std::endl;
    //     }
    
    // for(int b = (N+1)/2 + 3 ; b <= N-2 ; b +=3 )
    // // for(int b = (N+1)/2 + 1 ; b <= N-2 ; b +=3 )
    //     {
    //     psi.position(b); 

    //     auto l = leftLinkIndex(psi,b);
    //     auto s = siteIndex(psi,b);
    //     auto [U,S,V] = svd(psi(b),{l,s});
    //     auto u = commonIndex(U,S);

    //     Real SvN = 0.;
    //     for(auto n : range1(dim(u)))
    //         {
    //         auto Sn = elt(S,n,n);
    //         auto p = sqr(Sn);
    //         if(p > 1E-12) SvN += -p*log(p);
    //         }
    //     printfln("Across bond b=%d, SvN = %.10f",b,SvN);
    //     outfile7 << SvN << std::endl;
    //     }
    
    // outfile7.close();





   return 0 ;
    }