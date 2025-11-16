#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    MPO G1 ;
    MPO G2 ;
    int N = 180 ;
    int mid = N/2 ;
    auto sites = Fermion(N,{"ConserveQNs=", false});
    readFromFile("outputpsi/sites_N_120_sweep_20_t_0.030000_J_1_V_0.600000_D_0.018000",sites);
    MPS psi(sites);
    readFromFile<MPS>("outputpsi/N_120_sweep_20_t_0.030000_J_1_V_0.600000_D_0.018000", psi);


    // measuring majorana with poor efficiency , used only for testing
    // std::ofstream outfile("output");
    // for(int j = 3; j <= N; j +=3)
    //     {
    //     // //re-gauge psi to get ready to measure at position j
    //     auto ampo3 = AutoMPO(sites) ;
    //     ampo3 += 1, "C", 3, "C", j;
    //     G1 = toMPO(ampo3);
    //     auto ampo4 = AutoMPO(sites) ;
    //     ampo4 += -1, "C", 3, "Cdag", j;
    //     G2 = toMPO(ampo4);
    //     auto ampo5 = AutoMPO(sites) ;
    //     ampo5 += 1, "Cdag", 3, "C", j;
    //     G3 = toMPO(ampo5);
    //     auto ampo6 = AutoMPO(sites) ;
    //     ampo6 += -1, "Cdag", 3, "Cdag", j;
    //     G4 = toMPO(ampo6);
    //     auto wave1 = inner(psi , G1 , psi) ;
    //     auto wave2 = inner(psi , G2 , psi) ;
    //     auto wave3 = inner(psi , G3 , psi) ;
    //     auto wave4 = inner(psi , G4 , psi) ;
    //     auto majorana= wave1 + wave2 + wave3 + wave4 ;
    //     printfln("Site %d occupation: %.5f", j, majorana);
    //     outfile << majorana << std::endl;
    //     }
    //     //print(result) ; 
    // outfile.close();
    


    // std::ofstream outfile("/mnt/d/OneDrive/programs/vscodepython/marjorana/output3");
    // for(int j = 6; j <= N-4; j +=3)
    //     {
    //     // //re-gauge psi to get ready to measure at position j
    //     auto ampo3 = AutoMPO(sites) ;
    //     ampo3 += 1, "C", 3, "C", j;
    //     G1 = toMPO(ampo3);
    //     auto ampo4 = AutoMPO(sites) ;
    //     ampo4 += 1, "C", 3, "Cdag", j;
    //     G2 = toMPO(ampo4);
    //     auto wave1 = inner(psi , G1 , psi) ;
    //     auto wave2 = inner(psi , G2 , psi) ;
    //     auto majorana= wave1 + wave2  ;
    //     printfln("Site %d occupation: %.5f", j, majorana);
    //     outfile << majorana << std::endl;
    //     }
    //     //print(result) ; 
    // outfile.close();
     
    //calcualting left majorana <c_j gamma_3>
    // std::ofstream outfile3("/mnt/d/OneDrive/programs/vscodepython/marjorana/output3");
    // for(int j = 6; j < N-6; j+=3)       
    // {
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



    //calcualting  majorana <GAMMA_i gamma_2j+1>
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
        auto result = elt(A3Adagj) ; // - elt(A3Aj)  - elt(Adag3Adagj) + elt(Adag3Aj);
        outfile9  << result << std::endl;
        }
    outfile9.close();

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
    //     outfile4   << result << std::endl;
    //     }
    // outfile4.close();

    // // calculate entanglment entropy for orbital B
    // std::ofstream outfile5("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try1");
    // for(int b = 4 ; b < N-1; b +=3 )
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


    // // orbital C
    // std::ofstream outfile6("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try2");
    // for(int b = 3 ; b < N-1; b +=3 )
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

    // // orbital A
    // std::ofstream outfile7("/mnt/d/OneDrive/programs/vscodepython/marjorana/entanglement/try3");
    // for(int b = 2 ; b < N-1; b +=3 )
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