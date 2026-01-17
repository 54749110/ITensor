#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>
using namespace itensor;

int main(int argc, char* argv[])
    {
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");
    auto N = input.getInt("N");
    auto t = input.getReal("t");
    auto J = input.getReal("J");
    auto V = input.getReal("V");
    auto D = input.getReal("D");
    auto U = input.getReal("U");
    auto mu = input.getReal("mu");
    auto sites = Fermion(N,{"ConserveQNs=", false}); 
    auto E = input.getReal("E");
    auto Nsweep = input.getInt("nsweeps");

    readFromFile("outputpsi/sites_N_90_sweep_85_t_0.030000_J_1.000000_V_0.600000_D_0.018000_U_0.790000_mu_0.000000_E_0.000000",sites);
    MPS psi(sites);
    readFromFile<MPS>("outputpsi/N_90_sweep_85_t_0.030000_J_1.000000_V_0.600000_D_0.018000_U_0.790000_mu_0.000000_E_0.000000", psi);

    std::string filebase = "N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_U_" + std::to_string(U)
                        + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename_gamma3gammaj = std::string("outputgamma3gammaj/try2gamma3gammaj_") + filebase;
    std::ofstream outfile1(filename_gamma3gammaj, std::ios::out | std::ios::trunc);
    for(int j = 6; j < N-1; j+=3)
        {
        auto Adag3 = op(sites,"Adag",3);
        auto A3 = op(sites,"A",3);
        auto Aj = op(sites,"A",j);
        auto Adagj = op(sites,"Adag",j);

        // guage psi is a must for contracting left side
        psi.position(3) ;
        auto psidag = dag(psi);
        psidag.prime();
        auto li_1 = leftLinkIndex(psi,3);

        //constructing majorana operator using spinless fermion basis
        auto Adag3Aj = prime(psi(3),li_1)*Adag3*psidag(3);
        auto  A3Aj = prime(psi(3),li_1)*A3*psidag(3);
        auto Adag3Adagj = prime(psi(3),li_1)*Adag3*psidag(3);
        auto  A3Adagj = prime(psi(3),li_1)*A3*psidag(3);

        for(int k = 4; k < j; ++k)
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
        auto result = elt(Adag3Aj) - elt(A3Aj) - elt(Adag3Adagj) + elt(A3Adagj);
        outfile1  << result << std::endl;
        }
    outfile1.close();

   return 0 ;
    }