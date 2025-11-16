#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
{
    // Define the set of N values
    std::vector<int> N_values = {60, 66, 72, 90, 102, 120 , 132, 150 , 162, 180, 192, 210, 240, 270, 300, 360, 450};

    for(auto N : N_values)
    {
        MPO H;
        MPS psi0;
        MPO Hpbc;
        MPO n_j;
        MPO G1;
        MPO G2;
        MPO G3;
        MPO G4;

        // Initialize the site degrees of freedom
        auto sites = Fermion(N,{"ConserveQNs=", false}); 

        auto t= 0.03;
        auto J= 1;
        auto V= 0.1;
        auto mu = 0; 
        auto D = 0.6*t ;

        // Open boundary condition
        auto ampo = AutoMPO(sites);
        for (int j = 1; j <= N-3 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+3, "C", j;
            ampo += -t, "Cdag", j,   "C", j+3;
        }
        for(int j = 3; j <= N-2; j += 3) // electron hopping over orbitals
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

        for (int j = 1; j <= N-1 ; j += 3) //onsite V
        {
            ampo += V, "Cdag", j, "C", j;
        }

        for (int j = 2; j <= N ; j += 3) //onsite V
        {
            ampo += -V, "Cdag", j, "C", j;
        }

        ampo += -J, "Cdag", N, "C", N-2;
        ampo += -J, "Cdag", N-2,   "C", N;
        ampo += -J, "Cdag", N, "C", N-1;
        ampo += -J, "Cdag", N-1,   "C", N;
        
        // electron pairing
        for (int j = 1; j <= N-3 ; j += 1) //electron pairing
        {
            ampo += D, "Cdag", j+3, "Cdag", j;
            ampo += D, "C", j,   "C", j+3;
        }

        H = toMPO(ampo);

        auto state = InitState(sites);
        for(auto i : range1(N))
        {
            if(i%3 == 3) state.set(i,"Occ");
            else         state.set(i,"Emp");
        }
        psi0 = MPS(state);

        auto sweeps = Sweeps(20);
        sweeps.maxdim() = 10,20,40,100,200,400,800,1600,3200,6400;
        sweeps.cutoff() = 1E-10;
        sweeps.niter() = 2;
        sweeps.noise() = 1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,1E-9,1E-10,0.0;

        auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
        printfln("N = %d, Ground State Energy = %.10f", N, energy);

        // 马约拉纳波函数左右半支，不同格点数，外插用
    //     std::string filename = "/mnt/d/OneDrive/programs/vscodepython/marjorana/wave/output_N_" + std::to_string(N) + ".txt";
    //     std::ofstream outfile3(filename);
    //     for(int j = 6; j < N-1; j += 3)
    //     {
    //         auto Adag3 = op(sites,"Adag",3);
    //         auto A3 = op(sites,"A",3);
    //         auto Aj = op(sites,"A",j);

    //         psi.position(3);
    //         auto psidag = dag(psi);
    //         psidag.prime();
    //         auto li_1 = leftLinkIndex(psi,3);
    //         auto Adag3Aj = prime(psi(3),li_1)*Adag3*psidag(3);
    //         auto  A3Aj = prime(psi(3),li_1)*A3*psidag(3);
    //         for(int k = 4; k < j; ++k)
    //         {
    //             Adag3Aj *= psi(k);
    //             Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
    //             Adag3Aj *= psidag(k);
    //             A3Aj *= psi(k);
    //             A3Aj *= op(sites,"F",k); //Jordan-Wigner string
    //             A3Aj *= psidag(k);
    //         }
    //         auto lj = rightLinkIndex(psi,j);
    //         Adag3Aj  *= prime(psi(j),lj);
    //         Adag3Aj  *= Aj;
    //         Adag3Aj  *= psidag(j);
    //         A3Aj  *= prime(psi(j),lj);
    //         A3Aj  *= Aj;
    //         A3Aj  *= psidag(j);

    //         auto result = elt(Adag3Aj) - elt(A3Aj);
    //         outfile3 << result << std::endl;
    //     }
    //     outfile3.close();

        std::string filename = "/mnt/d/OneDrive/programs/vscodepython/marjorana/waveR/output_N_" + std::to_string(N) + ".txt";
        std::ofstream outfile4(filename);
        for(int j = 6; j < N-1; j +=3)
            {
            auto AdagN = op(sites,"Adag",N);
            auto AN = op(sites,"A",N);
            auto Aj = op(sites,"A",j);

            psi.position(j) ;
            auto psidag = dag(psi);
            psidag.prime();
            auto li_j = leftLinkIndex(psi,j);
            auto AjAdagN = prime(psi(j),li_j)*Aj*psidag(j);
            auto  AjAN = prime(psi(j),li_j)*Aj*psidag(j);
            for(int k = j+1 ; k < N; ++k)
                {
                AjAdagN *= psi(k);
                AjAdagN *= op(sites,"F",k); //Jordan-Wigner string
                AjAdagN *= psidag(k);
                AjAN *= psi(k);
                AjAN *= op(sites,"F",k); //Jordan-Wigner string
                AjAN *= psidag(k);
                }
            AjAdagN  *= psi(N) ;
            AjAdagN  *= AdagN;
            AjAdagN  *= psidag(N);
            AjAN  *= psi(N) ;
            AjAN  *= AN;
            AjAN  *= psidag(N);
        
            //consider JW-transformaton the second sign is plus ,first sign is minus
            auto result = - elt(AjAdagN) + elt(AjAN);
            outfile4  << result << std::endl;
            }
        outfile4.close();

    }

    return 0;
}