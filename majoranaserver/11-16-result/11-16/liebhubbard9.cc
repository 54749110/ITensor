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
    H = toMPO(ampo);

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



    // initial state
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
    std::string filename = std::string("run_outputs/outputpsi/") + filebase;             // main psi file
    std::string filename2 = std::string("run_outputs/outputpsi/sites_") + filebase;      // sites file
    std::string filename_energy = std::string("run_outputs/outputcheck/energy_") + filebase; // energy log
    std::string filename_ent = std::string("run_outputs/outputcheck/ent_") + filebase;     // entanglement log
    std::string filename_gamma3gammaj = std::string("run_outputs/outputgamma3gammaj/gamma3gammaj_") + filebase;
    std::string filename_gammajgammaN = std::string("run_outputs/outputgammajgammaN/gammajgammaN_") + filebase;
    std::string filename_gammaigammaj = std::string("run_outputs/outputgammaigammaj/gammaigammaj_") + filebase;
    std::string filename_cicj = std::string("run_outputs/outputcicj/cicj_") + filebase;
    std::string filename_cicdagj = std::string("run_outputs/outputcicdagj/cicdagj_") + filebase;
    std::string filename_ni = std::string("run_outputs/outputni/ni_") + filebase;
    std::string filename_ninj= std::string("run_outputs/outputninj/ninj_") + filebase;
    std::string filename_EE= std::string("run_outputs/outputEE/EE_") + filebase;
    std::string filename_dEdU= "run_outputs/outputdEdU/dEdU_N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename_densitycorr= "run_outputs/outputdensitycorr/densitycorr_" + filebase ;
    std::string filename_CDWorder = "run_outputs/outputCDWorder/CDWorder_N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V) + "_U_" + std::to_string(U)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu)  ;

                        
    // clear existing files
    {
    std::ofstream ofs2(filename2, std::ios::out | std::ios::trunc);
    }
    {
    std::ofstream ofs1(filename, std::ios::out | std::ios::trunc);
    }
    {
    std::ofstream ofs4(filename_dEdU, std::ios::out | std::ios::trunc);
    }
    {
    std::ofstream ofs3(filename_CDWorder, std::ios::out | std::ios::trunc);
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
            ampo4 += U, "Cdag", j, "C", j ,"Cdag", j+3, "C", j+3 ;
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
    if(N%6 == 1) 
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


    // monitor the wall-clock time of this program
    auto wall_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> wall_elapsed = wall_end - wall_start;
    double wall_seconds = wall_elapsed.count();
    printfln("Wall-clock time (s) = %.6f", wall_seconds);

    return 0 ;
}