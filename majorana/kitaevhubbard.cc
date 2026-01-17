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
            std::vector<Real> ent_line(Nsites,0.0);
            // make a non-const local copy of psi so we can call position
            MPS psi_local = psi();

            // typical EE calculation
            for(int b = 1; b <= Nsites-1; b += 1)
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
                ent_line[b-1] = SvN;
                }

            // append one line (one sweep) with Nsites numbers separated by spaces
            std::ofstream ofs(ent_outfile_, std::ios::out | std::ios::app);
            if(ofs)
                {
                for(int i = 0; i < Nsites; ++i)
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

    // Open boundary condition hamiltonian of kitaev-hubbard model
    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-1 ; j += 1) //electron hopping 
        {
            ampo += -1, "Cdag", j+1, "C", j;
            ampo += -1, "Cdag", j,   "C", j+1;
        }

    // electron pairing
    for (int j = 1; j <= N-1 ; j += 1) //electron pairing
        {
            ampo += 1, "Cdag", j+1, "Cdag", j;
            ampo += 1, "C", j,   "C", j+1;
    }

    for (int j = 1; j <= N ; j += 1) //chemical potential
        {
            ampo += mu, "Cdag", j, "C", j;
            ampo += -mu, "C", j, "Cdag", j;
        }

    for (int j = 1; j <= N-1 ; j += 1) //chemical potential
        {
            ampo += -U, "Cdag", j, "C", j;
            ampo += U, "C", j, "Cdag", j;
        }

    for (int j = 2; j <= N ; j += 1) //chemical potential
        {
            ampo += -U, "Cdag", j, "C", j;
            ampo += U, "C", j, "Cdag", j;
        }

    // interaction
    for (int j = 1; j <= N-1 ; j += 1) 
        {
            ampo += 4*U, "Cdag", j, "C", j ,"Cdag", j+1, "C", j+1 ;
    }

    H = toMPO(ampo);


    // initial state N must be 303 306
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%2 == 0) state.set(i,"Occ");
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
    std::string filebase = "KN_"
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
    std::string filename_dEdU= "outputdEdU/dEdU_KN_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename_GSenergy= "outputGSenergy/GZSenergy_KN_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D) + "_mu_" + std::to_string(mu) + "_E_" + std::to_string(E) ;
    std::string filename_densitycorr= "outputdensitycorr/densitycorr_" + filebase ;
    std::string filename_CDWorder = "outputCDWorder/CDWorder_KN_"
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

    // DMRG process , updating GS-energy and updating psi for each dmrg call 
    Args dmrg_args = Args("Quiet",true);
    MPS psi = psi0 ; 
    MyObserver myobs(psi,dmrg_args,filename_energy,filename_ent);
    auto energy = dmrg(psi,H,sweeps,myobs,dmrg_args);
    printfln("Initial energy = %.5f", inner(psi0,H,psi0) );
    printfln("\nGround State Energy = %.10f",energy);
   
    // store the final psi and siteset
    //writeToFile(filename2,sites);
    //writeToFile(filename, psi);

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

    // occupation number <n_j> as 3 aligned columns
    std::ofstream outfile6(filename_ni, std::ios::out | std::ios::trunc);
    // prepare three column index lists
    std::vector<int> col1; // 1,2,3,...
    for(int j = 1; j <= N; j += 1) col1.push_back(j);
    size_t maxRows3 = std::max({col1.size()});

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

        outfile6 << std::endl;
        }
    outfile6.close();
 
    // entanglement entropy 
    std::ofstream outfile8(filename_EE, std::ios::out | std::ios::trunc);
    for(int b = 1 ; b < N-1; b +=1 )
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
    for (int j = 1; j <= N-1 ; j += 1) //chemical potential
        {
            ampo4 += -1, "Cdag", j, "C", j;
            ampo4 += 1, "C", j, "Cdag", j;
        }

    for (int j = 2; j <= N ; j += 1) //chemical potential
        {
            ampo4 += -1, "Cdag", j, "C", j;
            ampo4 += 1, "C", j, "Cdag", j;
        }

    // interaction
    for (int j = 1; j <= N-1 ; j += 1) 
        {
            ampo4 += 4*1, "Cdag", j, "C", j ,"Cdag", j+1, "C", j+1 ;
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

    std::ofstream outfile12(filename_GSenergy, std::ios::out | std::ios::app);
    if(outfile12)
        {
        outfile12 << std::fixed << std::setprecision(12) << std::setw(22) << U
                 << ' ' << std::fixed << std::setprecision(12) << std::setw(22) << energy
                 << std::endl;
        }
    outfile12.close();

    // monitor the wall-clock time of this program
    auto wall_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> wall_elapsed = wall_end - wall_start;
    double wall_seconds = wall_elapsed.count();
    printfln("Wall-clock time (s) = %.6f", wall_seconds);

    return 0 ;
}