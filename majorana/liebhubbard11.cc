#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

using namespace itensor;

// the class defined below monitor the sweep process to justify convergence of DMRG by Energy and Entanglement entropy
// DMRGObserver class from ITensor to define MyObserver class
class MyObserver : public DMRGObserver
    {
    public:
    MyObserver(MPS const& psi, Args const& args = Args::global(), std::string energy_outfile = "sweep_energy.txt", std::string ent_outfile = "entanglement.txt")
        // syntax from Itensor to get information once every "dmrg" call is made.
        : DMRGObserver(psi,args), energy_outfile_(std::move(energy_outfile)), ent_outfile_(std::move(ent_outfile))
        {
        // clear files 
        std::ofstream ofs1(energy_outfile_, std::ios::out | std::ios::trunc);
        std::ofstream ofs2(ent_outfile_, std::ios::out | std::ios::trunc);

        }

    void measure(Args const& args) override
        {
        DMRGObserver::measure(args);
        // Only act at the end of a full sweep: HalfSweep==2 and AtBond==1 
        auto ha = args.getInt("HalfSweep",0);
        auto atbond = args.getInt("AtBond",0);
        int Nsites = length(psi());
        if(ha == 2 && atbond == 1)
            {
            // calling energy by Itensor syntax
            auto energy = args.getReal("Energy",0.);
            {
            std::ofstream ofs(energy_outfile_, std::ios::out | std::ios::app);
            if(ofs)
                {
                ofs << std::fixed << std::setprecision(12) << std::setw(22) << energy << std::endl;
                }
            }

            // compute entanglement entropy across internal bonds and emit a single line
            std::vector<Real> ent_line(Nsites/3,0.0);
            // make a non-const local copy of psi so we can call position
            MPS psi_local = psi();

            // typical EE calculation
            for(int b = 3; b <= Nsites-1; b += 3)
                {
                psi_local.position(b);
                auto l = leftLinkIndex(psi_local,b);
                auto s = siteIndex(psi_local,b);
                auto [U,S,V] = svd(psi_local(b),{l,s});
                auto u = commonIndex(U,S);
                Real SvN = 0.;
                for(auto n : range1(dim(u)))
                    {
                    auto Sn = elt(S,n,n);
                    auto p = sqr(Sn);
                    if(p > 1E-12) SvN += -p*log(p);
                    }
                ent_line[b/3-1] = SvN;
                }

            // append one line (one sweep) with Nsites numbers separated by spaces
            std::ofstream ofs(ent_outfile_, std::ios::out | std::ios::app);
            if(ofs)
                {
                for(int i = 0; i < Nsites/3; ++i)
                    {
                    if(i) ofs << ' ';
                    ofs << std::fixed << std::setprecision(12) << std::setw(22) << ent_line[i];
                    }
                ofs << std::endl;
                }
            }
        }

    private:
    std::string energy_outfile_;
    std::string ent_outfile_;
    };

int main(int argc, char* argv[])
{
    // start wall-clock timer
    auto wall_start = std::chrono::steady_clock::now();
    // read from input file
    // the command is ./liebhubbard2 liebinput
    println("//////////////////////////");
    println("Reading input file ......\n");
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    MPO H;
    MPS psi0;

    // attractive hubbard model
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    //  N mod 4 = 0 ; 
    auto N = input.getInt("N");
    auto t = input.getReal("t");
    auto J = input.getReal("J");
    auto V = input.getReal("V");
    auto D = input.getReal("D");
    auto U = input.getReal("U");
    auto mu = input.getReal("mu");
    auto sites = Fermion(N,{"ConserveQNs=", false}); 
    int mid = N/2 ; 
    auto E = input.getReal("E");

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

    for (int j = 1; j <= N ; j += 1) //chemical potential
        {
            ampo += -0.5*mu, "Cdag", j, "C", j;
           ampo += 0.5*mu, "C", j, "Cdag", j;
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

    // interaction
    for (int j = 1; j <= N-3 ; j += 1) 
        {
            ampo += U, "Cdag", j, "C", j ,"Cdag", j+3, "C", j+3 ;
    }

    // interaction2
    for (int j = 1; j <= N-3 ; j += 1) 
        {
            ampo += -U/2, "Cdag", j, "C", j ;
    }
    for (int j = 4; j <= N ; j += 1) 
        {
            ampo += -U/2, "Cdag", j, "C", j ;
    }


    // pinning field for odd N, to make CDW (1010……10101) , two negative pinning fields at both ends
    // pinning field for odd N, to make CDW (1010……1010) , negative and positive pinning field at each end
    if(N%2 == 0) 
        {
        ampo += -0.5*E, "Cdag", 1, "C", 1;
        ampo += 0.5*E, "C", 1, "Cdag", 1;
        ampo += -0.5*E, "Cdag", 2,"C", 2;
        ampo += 0.5*E, "C", 2, "Cdag", 2;
        ampo += -0.5*E, "Cdag", 3, "C", 3;
        ampo += 0.5*E, "C", 3, "Cdag", 3;
        ampo += 0.5*E, "Cdag", N, "C", N;
        ampo += -0.5*E, "C", N, "Cdag", N;
        ampo += 0.5*E, "Cdag", N-1, "C", N-1;
        ampo += -0.5*E, "C", N-1, "Cdag", N-1;
        ampo += 0.5*E, "Cdag", N-2, "C", N-2;
        ampo += -0.5*E, "C", N-2, "Cdag", N-2;
        }

    if(N%2 == 1) 
        {
        ampo += -0.5*E, "Cdag", 1, "C", 1;
        ampo += 0.5*E, "C", 1, "Cdag", 1;
        ampo += -0.5*E, "Cdag", 2,"C", 2;
        ampo += 0.5*E, "C", 2, "Cdag", 2;
        ampo += -0.5*E, "Cdag", 3, "C", 3;
        ampo += 0.5*E, "C", 3, "Cdag", 3;
        ampo += -0.5*E, "Cdag", N, "C", N;
        ampo += 0.5*E, "C", N, "Cdag", N;
        ampo += -0.5*E, "Cdag", N-1, "C", N-1;
        ampo += 0.5*E, "C", N-1, "Cdag", N-1;
        ampo += -0.5*E, "Cdag", N-2, "C", N-2;
        ampo += 0.5*E, "C", N-2, "Cdag", N-2;
        }

    H = toMPO(ampo);


    // initial state N must be 303 306
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%3 == 3) state.set(i,"Occ");
        else         state.set(i,"Emp");
        }
    psi0 = MPS(state);

    // sweep information
    auto Nsweep = input.getInt("nsweeps");
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(Nsweep,table);
    println(sweeps);


    // definition of output files
    // store both the MPS and MPS-siteset after DMRG for convenient calling
    std::string filebase = "N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_U_" + std::to_string(U)
                        + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename = std::string("outputpsi/") + filebase;             // main psi file
    std::string filename2 = std::string("outputpsi/sites_") + filebase;      // sites file
    std::string filename_energy = std::string("outputcheck/energy_") + filebase; // energy log
    std::string filename_ent = std::string("outputcheck/ent_") + filebase;     // entanglement log
    std::string filename_gamma3gammaj = std::string("outputgamma3gammaj/gamma3gammaj_") + filebase;
    std::string filename_gammajgammaN = std::string("outputgammajgammaN/gammajgammaN_") + filebase;
    std::string filename_gammaigammaj = std::string("outputgammaigammaj/gammaigammaj_") + filebase;
    std::string filename_cicj = std::string("outputcicj/cicj_") + filebase;
    std::string filename_cicdagj = std::string("outputcicdagj/cicdagj_") + filebase;
    std::string filename_ni = std::string("outputni/ni_") + filebase;
    std::string filename_ninj= std::string("outputninj/ninj_") + filebase;
    std::string filename_EE= std::string("outputEE/EE_") + filebase;
    std::string filename_pdw= std::string("outputpdw/pdw_") + filebase;
    std::string filename_dEdU= "outputdEdU/dEdU_N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename_GSenergy= "outputGSenergy/GZSenergy_N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename_densitycorr= "outputdensitycorr/densitycorr_" + filebase ;
    std::string filename_CDWorder = "outputCDWorder/CDWorder_N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V) + "_U_" + std::to_string(U)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu)  ;
    std::string filename_AdagB = std::string("output2cdagicj/AdagB_") + filebase;
    std::string filename_BdagA = std::string("output2cdagicj/BdagA_") + filebase;
    std::string filename_AdagC = std::string("output2cdagicj/AdagC_") + filebase;
    std::string filename_CdagA = std::string("output2cdagicj/CdagA_") + filebase;
    std::string filename_BdagC = std::string("output2cdagicj/BdagC_") + filebase;
    std::string filename_CdagB = std::string("output2cdagicj/CdagB_") + filebase;
    std::string filename_AdagA = std::string("output2cdagicj/AdagA_") + filebase;
    std::string filename_BdagB = std::string("output2cdagicj/BdagB_") + filebase;
    std::string filename_CdagC = std::string("output2cdagicj/CdagC_") + filebase;
    std::string filename_AA = std::string("output2cicj/AA_") + filebase;
    std::string filename_BB = std::string("output2cicj/BB_") + filebase;
    std::string filename_CC = std::string("output2cicj/CC_") + filebase;
    std::string filename_AB = std::string("output2cicj/AB_") + filebase;
    std::string filename_BA = std::string("output2cicj/BA_") + filebase;
    std::string filename_AC = std::string("output2cicj/AC_") + filebase;
    std::string filename_CA = std::string("output2cicj/CA_") + filebase;
    std::string filename_BC = std::string("output2cicj/BC_") + filebase;
    std::string filename_CB = std::string("output2cicj/CB_") + filebase;
    std::string filename_nAB = std::string("output2ninj/nAB_") + filebase;
    std::string filename_nAC = std::string("output2ninj/nAC_") + filebase;
    std::string filename_nBC = std::string("output2ninj/nBC_") + filebase;
    std::string filename_nAA = std::string("output2ninj/nAA_") + filebase;
    std::string filename_nBB = std::string("output2ninj/nBB_") + filebase;
    std::string filename_nCC = std::string("output2ninj/nCC_") + filebase;
    std::string filename_corrAB = std::string("output2corr/corrAB_") + filebase;
    std::string filename_corrAC = std::string("output2corr/corrAC_") + filebase;
    std::string filename_corrBC = std::string("output2corr/corrBC_") + filebase;
    std::string filename_corrAA = std::string("output2corr/corrAA_") + filebase;
    std::string filename_corrBB = std::string("output2corr/corrBB_") + filebase;
    std::string filename_corrCC = std::string("output2corr/corrCC_") + filebase;

                        
    // clear existing files
    {
    std::ofstream ofs2(filename2, std::ios::out | std::ios::trunc);
    }
    {
    std::ofstream ofs1(filename, std::ios::out | std::ios::trunc);
    }

    // DMRG process , updating GS-energy and updating psi for each dmrg call 
    Args dmrg_args = Args("Quiet",true);
    MPS psi = psi0 ; 
    MyObserver myobs(psi,dmrg_args,filename_energy,filename_ent);
    auto energy = dmrg(psi,H,sweeps,myobs,dmrg_args);
    printfln("Initial energy = %.5f", inner(psi0,H,psi0) );
    printfln("\nGround State Energy = %.10f",energy);
   
    // store the final psi and siteset
    writeToFile(filename2,sites);
    writeToFile(filename, psi);

    // read energy and entanglement files to check convergence
    auto check_convergence = [&](const std::string &energy_file, const std::string &ent_file)->bool
        {
        // convergence 10^-7
        const double tol = 1e-7;
        // energy
        std::vector<double> energies;
        {
        std::ifstream ifs(energy_file);
        if(!ifs) return false;
        double v;
        while(ifs >> v) energies.push_back(v);
        }
        if(energies.size() < 2) return false;
        double e_last = energies[energies.size()-1];
        double e_prev = energies[energies.size()-2];
        if(std::fabs(e_last - e_prev) >= tol) return false;

        // EE
        std::vector<std::vector<double>> ent_lines;
        {
        std::ifstream ifs(ent_file);
        if(!ifs) return false;
        std::string line;
        while(std::getline(ifs,line))
            {
            if(line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
            std::istringstream iss(line);
            std::vector<double> vals;
            double x;
            while(iss >> x) vals.push_back(x);
            if(!vals.empty()) ent_lines.push_back(std::move(vals));
            }
        }
        if(ent_lines.size() < 2) return false;
        const auto &last = ent_lines.back();
        const auto &prev = ent_lines[ent_lines.size()-2];
        if(last.size() != prev.size()) return false;
        for(size_t i = 0; i < last.size(); ++i)
            {
            if(std::fabs(last[i] - prev[i]) >= tol) return false;
            }
        return true;
        };

    if(check_convergence(filename_energy, filename_ent))
        {
        printfln("%s converge", filename.c_str());
        }
    else
        {
        printfln("%s notconverge", filename.c_str());
        }
    

    // the program below are the observations

    // calcualting left majorana <GAMMA_3 gamma_2J+1>
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

    //calcualting right majorana <c_j*i*gamma_N>
    std::ofstream outfile2(filename_gammajgammaN, std::ios::out | std::ios::trunc);
    for(int j = 6; j < N-1; j +=3)
        {
        auto AdagN = op(sites,"Adag",N);
        auto AN = op(sites,"A",N);
        auto Aj = op(sites,"A",j);
        auto Adagj = op(sites,"Adag",j);

        psi.position(j) ;
        auto psidag = dag(psi);
        psidag.prime();
        auto li_j = leftLinkIndex(psi,j);

        auto AjAdagN = prime(psi(j),li_j)*Aj*psidag(j);
        auto  AjAN = prime(psi(j),li_j)*Aj*psidag(j);
        auto AdagjAdagN = prime(psi(j),li_j)*Adagj*psidag(j);
        auto  AdagjAN = prime(psi(j),li_j)*Adagj*psidag(j);

        for(int k = j+1 ; k < N; ++k)
            {
            AjAdagN *= psi(k);
            AjAdagN *= op(sites,"F",k); //Jordan-Wigner string
            AjAdagN *= psidag(k);
            AjAN *= psi(k);
            AjAN *= op(sites,"F",k); //Jordan-Wigner string
            AjAN *= psidag(k);
            AdagjAdagN *= psi(k);
            AdagjAdagN *= op(sites,"F",k); //Jordan-Wigner string
            AdagjAdagN *= psidag(k);
            AdagjAN *= psi(k);
            AdagjAN *= op(sites,"F",k); //Jordan-Wigner string
            AdagjAN *= psidag(k);
            }

        AjAdagN  *= psi(N) ;
        AjAdagN  *= AdagN;
        AjAdagN  *= psidag(N);
        AjAN  *= psi(N) ;
        AjAN  *= AN;
        AjAN  *= psidag(N);
        AdagjAdagN  *= psi(N) ;
        AdagjAdagN  *= AdagN;
        AdagjAdagN  *= psidag(N);
        AdagjAN  *= psi(N) ;
        AdagjAN  *= AN;
        AdagjAN  *= psidag(N);
       
        //consider JW-transformaton the second sign is plus ,first sign is minus
        //origin : -cjcdagN + cjcN + cdagjcN - cdagjcdagN
        //now:    ajadagj - ajaN + adagjaN - adagjadagN 
        auto result =  elt(AjAdagN) - elt(AjAN) + elt(AdagjAN) - elt(AdagjAdagN) ;
        outfile2 << result << std::endl;
        }
    outfile2.close();

    // calcualting majorana matrix: columns are i = mid/2, mid/2+3, ..., mid
    // rows are j = mid/2, mid/2+3, ..., N-3. For j < i entries write "nan".
    std::ofstream outfile3(filename_gammaigammaj, std::ios::out | std::ios::trunc);
    std::ofstream outfile4(filename_cicj, std::ios::out | std::ios::trunc);
    std::ofstream outfile5(filename_cicdagj, std::ios::out | std::ios::trunc);
    int istart = mid/2;
    std::vector<int> ilist;
    for(int i = istart; i <= mid; i += 3) ilist.push_back(i);

    // determine number of rows by the longest column
    size_t maxRows = 0;
        for(auto i : ilist)
            {
            size_t rows = 0;
            // rows count for j = i+3, i+6, ..., <= N-3
            if(i + 3 <= N-3) rows = (N-3 - (i + 3))/3 + 1;
            if(rows > maxRows) maxRows = rows;
            }

    // for each row index r, compute j = ilist[c] + r*3 for each column c
        for(size_t r = 0; r < maxRows; ++r)
        {
        for(size_t c = 0; c < ilist.size(); ++c)
            {
            int i = ilist[c];
            if(c) outfile3 << ' ';
                // j starts from i+3
                int j = i + (static_cast<int>(r) + 1) * 3;
                if(j > N-3 || j <= i)
                {
                continue;
                }

            // compute observable for pair (i,j)
            auto Adagi = op(sites,"Adag",i);
            auto Ai = op(sites,"A",i);
            auto Aj = op(sites,"A",j);
            auto Adagj = op(sites,"Adag",j);

            psi.position(i);
            auto psidag = dag(psi);
            psidag.prime();
            auto li = leftLinkIndex(psi,i);

            auto AdagiAj = prime(psi(i),li)*Adagi*psidag(i);
            auto AiAj = prime(psi(i),li)*Ai*psidag(i);
            auto AdagiAdagj = prime(psi(i),li)*Adagi*psidag(i);
            auto AiAdagj = prime(psi(i),li)*Ai*psidag(i);

            for(int k = i+1; k < j; ++k)
                {
                AdagiAj *= psi(k);
                AdagiAj *= op(sites,"F",k);
                AdagiAj *= psidag(k);
                AiAj *= psi(k);
                AiAj *= op(sites,"F",k);
                AiAj *= psidag(k);
                AdagiAdagj *= psi(k);
                AdagiAdagj *= op(sites,"F",k);
                AdagiAdagj *= psidag(k);
                AiAdagj *= psi(k);
                AiAdagj *= op(sites,"F",k);
                AiAdagj *= psidag(k);
                }
            auto lj = rightLinkIndex(psi,j);

            AdagiAj  *= prime(psi(j),lj);
            AdagiAj  *= Aj;
            AdagiAj  *= psidag(j);
            AiAj  *= prime(psi(j),lj);
            AiAj  *= Aj;
            AiAj  *= psidag(j);
            AdagiAdagj  *= prime(psi(j),lj);
            AdagiAdagj  *= Adagj;
            AdagiAdagj  *= psidag(j);
            AiAdagj  *= prime(psi(j),lj);
            AiAdagj  *= Adagj;
            AiAdagj  *= psidag(j);

                auto result = elt(AdagiAj) - elt(AiAj) - elt(AdagiAdagj) + elt(AiAdagj);
                auto result2 = - elt(AiAj) ;
                auto result3 = - elt(AdagiAdagj);
                // align columns like MyObserver: a single space then setw(22)
                if(c)
                    {
                    outfile3 << ' ';
                    outfile4 << ' ';
                    outfile5 << ' ';
                    }
                outfile3 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                outfile4 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
                outfile5 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
            }
        outfile3 << std::endl;
        outfile4 << std::endl;
        outfile5 << std::endl;
        }
    outfile3.close();
    outfile4.close();
    outfile5.close();

    // occupation number <n_j> as 3 aligned columns
    std::ofstream outfile6(filename_ni, std::ios::out | std::ios::trunc);
    // prepare three column index lists
    std::vector<int> col1; // 3,6,9,...
    std::vector<int> col2; // 1,4,7,...
    std::vector<int> col3; // 2,5,8,...
    for(int j = 3; j <= N; j += 3) col1.push_back(j);
    for(int j = 1; j <= N; j += 3) col2.push_back(j);
    for(int j = 2; j <= N; j += 3) col3.push_back(j);
    size_t maxRows3 = std::max({col1.size(), col2.size(), col3.size()});

    // Nj1 will record the Nj values from the first column (column1) for each row
    std::vector<double> Nj1;

    for(size_t r = 0; r < maxRows3; ++r)
        {
        // column 1
        if(r < col1.size())
            {
            int j = col1[r];
            psi.position(j);
            auto ket = psi(j);
            auto bra = dag(prime(ket,"Site"));
            auto Njop = op(sites,"N",j);
            auto Nj = elt(bra*Njop*ket);
            outfile6 << std::fixed << std::setprecision(12) << std::setw(22) << Nj;
            Nj1.push_back(Nj);
            }
        else outfile6 << std::setw(22) << "nan";

        // column 2
        if(r < col2.size())
            {
            int j = col2[r];
            psi.position(j);
            auto ket = psi(j);
            auto bra = dag(prime(ket,"Site"));
            auto Njop = op(sites,"N",j);
            auto Nj = elt(bra*Njop*ket);
            outfile6 << ' ' << std::fixed << std::setprecision(12) << std::setw(22) << Nj;
            }
        else outfile6 << ' ' << std::setw(22) << "nan";

        // column 3
        if(r < col3.size())
            {
            int j = col3[r];
            psi.position(j);
            auto ket = psi(j);
            auto bra = dag(prime(ket,"Site"));
            auto Njop = op(sites,"N",j);
            auto Nj = elt(bra*Njop*ket);
            outfile6 << ' ' << std::fixed << std::setprecision(12) << std::setw(22) << Nj;
            }
        else outfile6 << ' ' << std::setw(22) << "nan";

        outfile6 << std::endl;
        }
    outfile6.close();

    MPO G1 ;
    // using autoMPO to calculate ninj
    std::ofstream outfile7(filename_ninj, std::ios::out | std::ios::trunc);
    int istart2 = 3;
    std::vector<int> ilist2;
    for(int i = istart2; i <= mid; i += 3) ilist2.push_back(i);

    // determine number of rows by the longest column
    size_t maxRows2 = 0;
        for(auto i : ilist2)
            {
            size_t rows2 = 0;
            // rows count for j = i+3, i+6, ..., <= N-3
            if(i + 3 <= N-3) rows2 = (N-3 - (i + 3))/3 + 1;
            if(rows2 > maxRows2) maxRows2 = rows2;
            }

    // for each row index r, compute j = ilist[c] + r*3 for each column c
        for(size_t r = 0; r < maxRows2; ++r)
        {
        for(size_t c = 0; c < ilist2.size(); ++c)
            {
            int i = ilist2[c];
            if(c) outfile7 << ' ';
                // j starts from i+3
                int j = i + (static_cast<int>(r) + 1) * 3;
                if(j > N-3 || j <= i)
                {
                continue;
                }

                // using MPO to represent ninj
                auto ampo3 = AutoMPO(sites) ;
                ampo3 += 1, "Cdag", i, "C", i , "Cdag", j, "C", j ;
                // build MPO from AutoMPO so its length matches the MPS
                G1 = toMPO(ampo3);
                auto wave1 = inner(psi , G1 , psi) ;
            if(c)
                    {
                    outfile7 << ' ';
                    }
                outfile7 << std::fixed << std::setprecision(12) << std::setw(22) << wave1;
            }
        outfile7 << std::endl;
        }
    outfile7.close();
 
    // entanglement entropy for orbital B
    std::ofstream outfile8(filename_EE, std::ios::out | std::ios::trunc);
    for(int b = 3 ; b < N-1; b +=3 )
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
        outfile8 << SvN << std::endl;
        }
    outfile8.close();

    // Feynman-Hellman theorem , partial|psi(U0)> / partialU = 0
    // phase transition partialEgs/partialU |(U0) = <psi(U0)|dH/dU|psi(U0)> = <psi(U0)| sum_i n_i n_{i+3} |psi(U0)>
    // calculate via MPO sum_i n_i n_{i+3} over G.S. wave function
    MPO G2 ;
    // ATTENTION ! the former file_name definition record U-wave2 for the same mu 
    std::ofstream outfile9(filename_dEdU, std::ios::out | std::ios::app);
    auto ampo4 = AutoMPO(sites) ;
    for (int j = 1; j <= N-3 ; j += 1) 
        {
            ampo4 += 1, "Cdag", j, "C", j ,"Cdag", j+3, "C", j+3 ;
    }
    G2 = toMPO(ampo4);
    auto wave2 = inner(psi , G2 , psi) ;
    if(outfile9)
        {
        outfile9 << std::fixed << std::setprecision(12) << std::setw(22) << U
                 << ' ' << std::fixed << std::setprecision(12) << std::setw(22) << wave2
                 << std::endl;
        }
    outfile9.close();

    // calcualting density correlation function <n_i n_j> - <n_i><n_j>
    // prepare <ni><nj> for different r = j - i
    // N=300, From averaging (25,26)……(50,51) to (25,50)……(50,75) as Nj2 elements length 25
    // N=303, From averaging (25,26)……(50,51) to (25,50)……(50,75) as Nj2 elements lenght 25
    // N=306, From averaging (25,26)……(51,52) to (25,51)……(51,77) as Nj2 elements length 26
    // N=309, From averaging (25,26)……(51,52) to (25,51)……(51,77) as Nj2 elements length 26
    std::vector<double> Nj2;
    int M = static_cast<int>(N/6); // floor
    int mid2 = static_cast<int>(N/12);
    for(int r = 1; r <= M-mid2; ++r)
        {
        double s = 0.0;
        for(int p = mid2; p <= M; ++p)
            {
            s += Nj1[p-1] * Nj1[p-1 + r];
            }
        Nj2.push_back(s/(M - mid2));
        }
 
    // prepare <ninj> for different r = j - i
    // N=300, From averaging (75,78)……(150,153) to (75,150)……(150,225) as Nj2 elements length 25
    // N=303, From averaging (75,78)……(150,153) to (75,150)……(150,225) as Nj2 elements lenght 25
    // N=306, From averaging (75,78)……(153,156) to (75,153)……(153,231) as Nj2 elements length 26
    // N=309, From averaging (75,78)……(153,156) to (75,153)……(153,231) as Nj2 elements length 26
    std::vector<double> ninj;
    int imax = 3*static_cast<int>(N/6);
    int imin = 3*static_cast<int>(N/12);
    for(int r = 3; r <= imax-imin; r += 3)
        {
        double s = 0.0;
        for(int p = imin; p <= imax; p += 3)
                {
                auto ampo_tmp = AutoMPO(sites);
                ampo_tmp += 1, "Cdag", p, "C", p, "Cdag", p+r, "C", p+r;
                auto Gtmp = toMPO(ampo_tmp);
                auto w = inner(psi, Gtmp, psi);
                s += w;
                }
            ninj.push_back(3*s/(imax - imin));
        }
    

    // calcualting density correlation function <n_i n_j> - <n_i><n_j>
    std::ofstream outfile10(filename_densitycorr, std::ios::out | std::ios::trunc);
    std::vector<double> density;
    for(size_t i = 0; i < Nj2.size(); ++i)
        {
        density.push_back(ninj[i] - Nj2[i]);
        }
    for(const auto& val : density)
        {
        outfile10 << std::fixed << std::setprecision(12) << std::setw(22) << val << std::endl;
        }
    outfile10.close();


   // calculating order parameter O = 1/N * sum((-1)^i*n_i) for different pinning field E
    std::ofstream outfile11(filename_CDWorder, std::ios::out | std::ios::app);
    double orderni = 0.0;
    if(N%6 == 0) 
        {
        for(int p = 1; p <= N/3; ++p)
            {
            orderni += 3*Nj1[p-1]*pow(-1,p)/N; 
            }
        }
    if(N%6 == 3) 
        {
        for(int p = 1; p <= N/3-1; ++p)
            {
            orderni += 3*Nj1[p-1]*pow(-1,p)/N; 
            }
        }
    if(outfile11)
        {
        outfile11 << std::fixed << std::setprecision(12) << std::setw(22) << E
                 << ' ' << std::fixed << std::setprecision(12) << std::setw(22) << orderni
                 << std::endl;
        }
    outfile11.close();

    std::ofstream outfile12(filename_GSenergy, std::ios::out | std::ios::app);
    if(outfile12)
        {
        outfile12 << std::fixed << std::setprecision(12) << std::setw(22) << U
                 << ' ' << std::fixed << std::setprecision(12) << std::setw(22) << energy
                 << std::endl;
        }
    outfile12.close();


    // calculating cdagAicBj to explore whether there is symmetry breaking between orbitals 
    // calculate cicj for any i, j , then classify them into different files according to orbitals
    // orbital A with V : i mod 3 = 1 ;orbital C with -V : i mod 3 = 2 ;orbital B with 0 : i mod 3 = 0
    // AB and BA are different bonds left and right to B
    std::ofstream outfile13(filename_AdagB, std::ios::out | std::ios::trunc);
    std::ofstream outfile37(filename_BdagA, std::ios::out | std::ios::trunc);
    std::ofstream outfile14(filename_AdagC, std::ios::out | std::ios::trunc);
    std::ofstream outfile38(filename_CdagA, std::ios::out | std::ios::trunc);
    std::ofstream outfile15(filename_BdagC, std::ios::out | std::ios::trunc);
    std::ofstream outfile39(filename_CdagB, std::ios::out | std::ios::trunc);
    std::ofstream outfile16(filename_AdagA, std::ios::out | std::ios::trunc);
    std::ofstream outfile17(filename_BdagB, std::ios::out | std::ios::trunc);
    std::ofstream outfile18(filename_CdagC, std::ios::out | std::ios::trunc);
    std::ofstream outfile19(filename_AB, std::ios::out | std::ios::trunc);
    std::ofstream outfile40(filename_BA, std::ios::out | std::ios::trunc);
    std::ofstream outfile20(filename_AC, std::ios::out | std::ios::trunc);
    std::ofstream outfile41(filename_CA, std::ios::out | std::ios::trunc);
    std::ofstream outfile21(filename_BC, std::ios::out | std::ios::trunc);
    std::ofstream outfile42(filename_CB, std::ios::out | std::ios::trunc);
    std::ofstream outfile22(filename_AA, std::ios::out | std::ios::trunc);
    std::ofstream outfile23(filename_BB, std::ios::out | std::ios::trunc);
    std::ofstream outfile24(filename_CC, std::ios::out | std::ios::trunc);

    for(int i = 4; i <= N-4; i += 1)
    {
        // 使用标志位记录每个文件在当前i循环中是否接收了数据
        bool hasData13 = false, hasData14 = false, hasData15 = false;
        bool hasData16 = false, hasData17 = false, hasData18 = false;
        bool hasData19 = false, hasData20 = false, hasData21 = false;
        bool hasData22 = false, hasData23 = false, hasData24 = false;
        bool hasData37 = false, hasData38 = false, hasData39 = false;
        bool hasData40 = false, hasData41 = false, hasData42 = false;
        
        // 遍历j之前先不输出空格
        for(int j = i+1; j <= N-3; j += 1)
        {
            // compute observable for pair (i,j)
            auto Adagi = op(sites,"Adag",i);
            auto Ai = op(sites,"A",i);
            auto Aj = op(sites,"A",j);
            auto Adagj = op(sites,"Adag",j);

            psi.position(i);
            auto psidag = dag(psi);
            psidag.prime();
            auto li = leftLinkIndex(psi,i);

            auto AdagiAj = prime(psi(i),li)*Adagi*psidag(i);
            auto AiAj = prime(psi(i),li)*Ai*psidag(i);
            auto AdagiAdagj = prime(psi(i),li)*Adagi*psidag(i);
            auto AiAdagj = prime(psi(i),li)*Ai*psidag(i);

            if (j - i > 1)
            {
                for(int k = i+1; k < j; ++k)
                {
                    AdagiAj *= psi(k);
                    AdagiAj *= op(sites,"F",k);
                    AdagiAj *= psidag(k);
                    AiAj *= psi(k);
                    AiAj *= op(sites,"F",k);
                    AiAj *= psidag(k);
                    AdagiAdagj *= psi(k);
                    AdagiAdagj *= op(sites,"F",k);
                    AdagiAdagj *= psidag(k);
                    AiAdagj *= psi(k);
                    AiAdagj *= op(sites,"F",k);
                    AiAdagj *= psidag(k);
                }
            }
            
            auto lj = rightLinkIndex(psi,j);

            AdagiAj  *= prime(psi(j),lj);
            AdagiAj  *= Aj;
            AdagiAj  *= psidag(j);
            AiAj  *= prime(psi(j),lj);
            AiAj  *= Aj;
            AiAj  *= psidag(j);
            AdagiAdagj  *= prime(psi(j),lj);
            AdagiAdagj  *= Adagj;
            AdagiAdagj  *= psidag(j);
            AiAdagj  *= prime(psi(j),lj);
            AiAdagj  *= Adagj;
            AiAdagj  *= psidag(j);

            auto result = elt(AdagiAj);
            auto result2 = - elt(AiAj);
            
            // 根据i,j的模值输出到对应的文件
            if (i % 3 == 1 && j % 3 == 0) // AB 
            {
                // 如果是第一个数据，先输出空格
                if (!hasData13) 
                {
                    outfile13 << ' ';
                    hasData13 = true;
                }
                outfile13 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData19) 
                {
                    outfile19 << ' ';
                    hasData19 = true;
                }
                outfile19 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 1 && j % 3 == 2) // AC 
            {
                if (!hasData14) 
                {
                    outfile14 << ' ';
                    hasData14 = true;
                }
                outfile14 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData20) 
                {
                    outfile20 << ' ';
                    hasData20 = true;
                }
                outfile20 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 0 && j % 3 == 2) // BC 
            {
                if (!hasData15) 
                {
                    outfile15 << ' ';
                    hasData15 = true;
                }
                outfile15 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData21) 
                {
                    outfile21 << ' ';
                    hasData21 = true;
                }
                outfile21 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 1 && j % 3 == 1) // AA 
            {
                if (!hasData16) 
                {
                    outfile16 << ' ';
                    hasData16 = true;
                }
                outfile16 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData22) 
                {
                    outfile22 << ' ';
                    hasData22 = true;
                }
                outfile22 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 0 && j % 3 == 0) // BB 
            {
                if (!hasData17) 
                {
                    outfile17 << ' ';
                    hasData17 = true;
                }
                outfile17 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData23) 
                {
                    outfile23 << ' ';
                    hasData23 = true;
                }
                outfile23 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 2 && j % 3 == 2) // CC
            {
                if (!hasData18) 
                {
                    outfile18 << ' ';
                    hasData18 = true;
                }
                outfile18 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData24) 
                {
                    outfile24 << ' ';
                    hasData24 = true;
                }
                outfile24 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 0 && j % 3 == 1) // BA
            {
                if (!hasData37) 
                {
                    outfile37 << ' ';
                    hasData37 = true;
                }
                outfile37 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData40) 
                {
                    outfile40 << ' ';
                    hasData40 = true;
                }
                outfile40 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 2 && j % 3 == 1) // CA
            {
                if (!hasData38) 
                {
                    outfile38 << ' ';
                    hasData38 = true;
                }
                outfile38 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData41) 
                {
                    outfile41 << ' ';
                    hasData41 = true;
                }
                outfile41 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            else if (i % 3 == 2 && j % 3 == 0) // CB
            {
                if (!hasData39) 
                {
                    outfile39 << ' ';
                    hasData39 = true;
                }
                outfile39 << std::fixed << std::setprecision(12) << std::setw(22) << result;
                
                if (!hasData42) 
                {
                    outfile42 << ' ';
                    hasData42 = true;
                }
                outfile42 << std::fixed << std::setprecision(12) << std::setw(22) << result2;
            }
            
        }
        
        // 只有接收了数据的文件才换行
        if (hasData13) outfile13 << std::endl;
        if (hasData14) outfile14 << std::endl;
        if (hasData15) outfile15 << std::endl;
        if (hasData16) outfile16 << std::endl;
        if (hasData17) outfile17 << std::endl;
        if (hasData18) outfile18 << std::endl;
        if (hasData19) outfile19 << std::endl;
        if (hasData20) outfile20 << std::endl;
        if (hasData21) outfile21 << std::endl;
        if (hasData22) outfile22 << std::endl;
        if (hasData23) outfile23 << std::endl;
        if (hasData24) outfile24 << std::endl;
        if (hasData37) outfile37 << std::endl;
        if (hasData38) outfile38 << std::endl;
        if (hasData39) outfile39 << std::endl;
        if (hasData40) outfile40 << std::endl;
        if (hasData41) outfile41 << std::endl;
        if (hasData42) outfile42 << std::endl;
    }

    outfile13.close();
    outfile14.close();
    outfile15.close();
    outfile16.close();
    outfile17.close();
    outfile18.close();
    outfile19.close();
    outfile20.close();
    outfile21.close();
    outfile22.close();
    outfile23.close();
    outfile24.close();
    outfile37.close();
    outfile38.close();
    outfile39.close();
    outfile40.close();
    outfile41.close();
    outfile42.close();

    // calaulating nAinBj and density correlation to explore whether there is symmetry breaking between orbitals 
    // calculate ninj and ni for any i, j , then classify them into different files according to orbitals
    // orbital A with V : i mod 3 = 1 ;orbital C with -V : i mod 3 = 2 ;orbital B with 0 : i mod 3 = 0  
    std::ofstream outfile25(filename_nAB, std::ios::out | std::ios::trunc);
    std::ofstream outfile26(filename_nAC, std::ios::out | std::ios::trunc);
    std::ofstream outfile27(filename_nBC, std::ios::out | std::ios::trunc);
    std::ofstream outfile28(filename_nAA, std::ios::out | std::ios::trunc);
    std::ofstream outfile29(filename_nBB, std::ios::out | std::ios::trunc);
    std::ofstream outfile30(filename_nCC, std::ios::out | std::ios::trunc);
    std::ofstream outfile31(filename_corrAB, std::ios::out | std::ios::trunc);
    std::ofstream outfile32(filename_corrAC, std::ios::out | std::ios::trunc);
    std::ofstream outfile33(filename_corrBC, std::ios::out | std::ios::trunc);
    std::ofstream outfile34(filename_corrAA, std::ios::out | std::ios::trunc);
    std::ofstream outfile35(filename_corrBB, std::ios::out | std::ios::trunc);
    std::ofstream outfile36(filename_corrCC, std::ios::out | std::ios::trunc);
    std::vector<int> occupy;
    for(int i = N/4; i <= N/2; i += 1)
        {
        psi.position(i);
        auto ket = psi(i);
        auto bra = dag(prime(ket,"Site"));
        auto Njop = op(sites,"N",i);
        auto Nj = elt(bra*Njop*ket);
        occupy.push_back(Nj);
        }
    for(int i = N/4; i <= N/2; i += 1)
        {

        bool hasData25 = false, hasData26 = false, hasData27 = false;
        bool hasData28 = false, hasData29 = false, hasData30 = false;
        bool hasData31 = false, hasData32 = false, hasData33 = false;
        bool hasData34 = false, hasData35 = false, hasData36 = false;

        for(int j = i; j <= N/2; j += 1)
            {
                auto ampo_tmp = AutoMPO(sites);
                ampo_tmp += 1, "Cdag", i, "C", i, "Cdag", j, "C", j;
                auto Gtmp = toMPO(ampo_tmp);
                auto result3 = inner(psi, Gtmp, psi);
                auto density_corr2 = result3 - occupy[i-1]*occupy[j-1];
                if (i % 3 == 1 && j % 3 == 0) // AB 
                {
                    // 如果是第一个数据，先输出空格
                    if (!hasData25) 
                    {
                        outfile25 << ' ';
                        hasData25 = true;
                    }
                    outfile25 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
                    
                    if (!hasData31) 
                    {
                        outfile31 << ' ';
                        hasData31 = true;
                    }
                    outfile31 << std::fixed << std::setprecision(12) << std::setw(22) << density_corr2;
                }
                else if (i % 3 == 1 && j % 3 == 2) // AC 
                {
                    if (!hasData26) 
                    {
                        outfile26 << ' ';
                        hasData26 = true;
                    }
                    outfile26 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
                    
                    if (!hasData32) 
                    {
                        outfile32 << ' ';
                        hasData32 = true;
                    }
                    outfile32 << std::fixed << std::setprecision(12) << std::setw(22) << density_corr2;
                }
                else if (i % 3 == 0 && j % 3 == 2) // BC 
                {
                    if (!hasData27) 
                    {
                        outfile27 << ' ';
                        hasData27 = true;
                    }
                    outfile27 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
                    
                    if (!hasData33) 
                    {
                        outfile33 << ' ';
                        hasData33 = true;
                    }
                    outfile33 << std::fixed << std::setprecision(12) << std::setw(22) << density_corr2;
                }
                else if (i % 3 == 1 && j % 3 == 1) // AA 
                {
                    if (!hasData28) 
                    {
                        outfile28 << ' ';
                        hasData28 = true;
                    }
                    outfile28 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
                    
                    if (!hasData34) 
                    {
                        outfile34 << ' ';
                        hasData34 = true;
                    }
                    outfile34 << std::fixed << std::setprecision(12) << std::setw(22) << density_corr2;
                }
                else if (i % 3 == 0 && j % 3 == 0) // BB 
                {
                    if (!hasData29) 
                    {
                        outfile29 << ' ';
                        hasData29 = true;
                    }
                    outfile29 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
                    
                    if (!hasData35) 
                    {
                        outfile35 << ' ';
                        hasData35 = true;
                    }
                    outfile35 << std::fixed << std::setprecision(12) << std::setw(22) << density_corr2;
                }
                else if (i % 3 == 2 && j % 3 == 2) // CC 
                {
                    if (!hasData30) 
                    {
                        outfile30 << ' ';
                        hasData30 = true;
                    }
                    outfile30 << std::fixed << std::setprecision(12) << std::setw(22) << result3;
                    
                    if (!hasData36) 
                    {
                        outfile36 << ' ';
                        hasData36 = true;
                    }
                    outfile36 << std::fixed << std::setprecision(12) << std::setw(22) << density_corr2;
                }               
            }

            // 只有接收了数据的文件才换行
            if (hasData25) outfile25 << std::endl;
            if (hasData26) outfile26 << std::endl;
            if (hasData27) outfile27 << std::endl;
            if (hasData28) outfile28 << std::endl;
            if (hasData29) outfile29 << std::endl;
            if (hasData30) outfile30 << std::endl;
            if (hasData31) outfile31 << std::endl;
            if (hasData32) outfile32 << std::endl;
            if (hasData33) outfile33 << std::endl;
            if (hasData34) outfile34 << std::endl;
            if (hasData35) outfile35 << std::endl;
            if (hasData36) outfile36 << std::endl;
        }

        outfile25.close();
        outfile26.close();
        outfile27.close();
        outfile28.close();
        outfile29.close();
        outfile30.close();
        outfile31.close();
        outfile32.close();
        outfile33.close();
        outfile34.close();
        outfile35.close();
        outfile36.close();


    // monitor the wall-clock time of this program
    auto wall_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> wall_elapsed = wall_end - wall_start;
    double wall_seconds = wall_elapsed.count();
    printfln("Wall-clock time (s) = %.6f", wall_seconds);

    return 0 ;
}