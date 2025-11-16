#include "itensor/all.h"
#include "itensor/util/print_macro.h"
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
    // read from input file
    // the command is ./liebhubbard2 liebinput
    println("//////////////////////////");
    println("Reading input file ......\n");
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    MPO H;
    MPS psi0;
    MPO Hpbc ;
    MPO n_j ;
    MPO G1;
    MPO G2;
    MPO G3;
    MPO G4;

    // attractive hubbard model
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    auto N = input.getInt("N");
    auto t = input.getReal("t");
    auto J = input.getReal("J");
    auto V = input.getReal("V");
    auto D = input.getReal("D");
    auto sites = Fermion(N,{"ConserveQNs=", false}); 

    // auto t= 0.03;
    // auto J= 1;
    // auto V= 0.6;
    // auto mu = 0; 
    // auto D = 0.6*t ;
    //int Nsweep = 20 ;
    int mid = N/4 ;

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
        
    for (int j = 1; j <= N-1 ; j += 3) //onsite V 这里几项我尝试配做particle-hole对称的形式，去掉了一个常数项
        {
            ampo += V, "Cdag", j, "C", j;
    }

    for (int j = 2; j <= N ; j += 3) //onsite V
        {
            ampo += -V, "Cdag", j, "C", j;
    }

    //for (int j = 1; j <= N ; j += 1) //chemical potential
    //    {
    //        ampo += -0.5*mu, "Cdag", j, "C", j;
    //       ampo += 0.5*mu, "C", j, "Cdag", j;
    //}

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

    
    // sweep
    // auto sweeps = Sweeps(Nsweep);
    // sweeps.maxdim() = 10,20,40,100,200,400,800,1600,3200,6400;
    // sweeps.cutoff() = 1E-10;
    // sweeps.niter() = 2;
    // sweeps.noise() = 1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,1E-9,1E-10,0.0;

    // store the monitoring energy and entanglement for each sweep
    // store both the MPS and MPS-siteset after DMRG for convenient calling
    // Use a single file base and derive related filenames from it
    std::string filebase = "N_"
                        + std::to_string(N) + "_sweep_" + std::to_string(Nsweep) + "_t_" + std::to_string(t)
                        + "_J_" + std::to_string(J) + "_V_" + std::to_string(V)
                        + "_D_" + std::to_string(D);
    std::string filename = std::string("outputpsi/") + filebase;             // main psi file
    std::string filename2 = std::string("outputpsi/sites_") + filebase;      // sites file
    std::string filename_energy = std::string("outputcheck/energy_") + filebase; // energy log
    std::string filename_ent = std::string("outputcheck/ent_") + filebase;     // entanglement log
    // store both the siteindex and psi

                        
    // clear existing files
    {
    std::ofstream ofs2(filename2, std::ios::out | std::ios::trunc);
    }
    {
    std::ofstream ofs1(filename, std::ios::out | std::ios::trunc);
    }
    

    // DMRG process , updating calculate energy and update psi for each dmrg call 
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
    

    // //measuring majorana with poor efficiency , used only for testing
    // std::ofstream outfile("/mnt/d/OneDrive/programs/vscodepython/marjorana/output6");
    // for(int j = 6; j <= N; j +=1)
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
    

    //calcualting left majorana <GAMMA_3 gamma_2J+1>
    // std::ofstream outfile8("/mnt/d/OneDrive/programs/vscodepython/marjorana/output5");
    // for(int j = 6; j < N-1; j+=3)
    //     {
    //     auto Adag3 = op(sites,"Adag",3);
    //     auto A3 = op(sites,"A",3);
    //     auto Aj = op(sites,"A",j);
    //     auto Adagj = op(sites,"Adag",j);


    //     // guage psi is a must for contracting left side
    //     psi.position(3) ;
    //     auto psidag = dag(psi);
    //     psidag.prime();
    //     auto li_1 = leftLinkIndex(psi,3);

    //     //constructing majorana operator using spinless fermion basis
    //     auto Adag3Aj = prime(psi(3),li_1)*Adag3*psidag(3);
    //     auto  A3Aj = prime(psi(3),li_1)*A3*psidag(3);
    //     auto Adag3Adagj = prime(psi(3),li_1)*Adag3*psidag(3);
    //     auto  A3Adagj = prime(psi(3),li_1)*A3*psidag(3);

    //     for(int k = 4; k < j; ++k)
    //         {
    //         Adag3Aj *= psi(k);
    //         Adag3Aj *= op(sites,"F",k); //Jordan-Wigner string
    //         Adag3Aj *= psidag(k);
    //         A3Aj *= psi(k);
    //         A3Aj *= op(sites,"F",k); //Jordan-Wigner string
    //         A3Aj *= psidag(k);
    //         Adag3Adagj *= psi(k);
    //         Adag3Adagj *= op(sites,"F",k); //Jordan-Wigner string
    //         Adag3Adagj *= psidag(k);
    //         A3Adagj *= psi(k);
    //         A3Adagj *= op(sites,"F",k); //Jordan-Wigner string
    //         A3Adagj *= psidag(k);
    //         }
    //     auto lj = rightLinkIndex(psi,j);

    //     Adag3Aj  *= prime(psi(j),lj);
    //     Adag3Aj  *= Aj;
    //     Adag3Aj  *= psidag(j);
    //     A3Aj  *= prime(psi(j),lj);
    //     A3Aj  *= Aj;
    //     A3Aj  *= psidag(j);
    //     Adag3Adagj  *= prime(psi(j),lj);
    //     Adag3Adagj  *= Adagj;
    //     Adag3Adagj  *= psidag(j);
    //     A3Adagj  *= prime(psi(j),lj);
    //     A3Adagj  *= Adagj;
    //     A3Adagj  *= psidag(j);

    //     //consider JW-transformaton 
    //     //origin : c3cj - c3cdagj + cdag3cj - cdag3cdagj
    //     //now:    -a3aj + a3adagj + adag3aj - adag3agdaj 
    //     auto result = elt(Adag3Aj) - elt(A3Aj) - elt(Adag3Adagj) + elt(A3Adagj);
    //     outfile8  << result << std::endl;
    //     }
    // outfile8.close();
    

    //calcualting right majorana <c_j*i*gamma_N>
    // std::ofstream outfile10("/mnt/d/OneDrive/programs/vscodepython/marjorana/output9");
    // for(int j = 6; j < N-1; j +=3)
    //     {
    //     auto AdagN = op(sites,"Adag",N);
    //     auto AN = op(sites,"A",N);
    //     auto Aj = op(sites,"A",j);
    //     auto Adagj = op(sites,"Adag",j);

    //     psi.position(j) ;
    //     auto psidag = dag(psi);
    //     psidag.prime();
    //     auto li_j = leftLinkIndex(psi,j);

    //     auto AjAdagN = prime(psi(j),li_j)*Aj*psidag(j);
    //     auto  AjAN = prime(psi(j),li_j)*Aj*psidag(j);
    //     auto AdagjAdagN = prime(psi(j),li_j)*Adagj*psidag(j);
    //     auto  AdagjAN = prime(psi(j),li_j)*Adagj*psidag(j);

    //     for(int k = j+1 ; k < N; ++k)
    //         {
    //         AjAdagN *= psi(k);
    //         AjAdagN *= op(sites,"F",k); //Jordan-Wigner string
    //         AjAdagN *= psidag(k);
    //         AjAN *= psi(k);
    //         AjAN *= op(sites,"F",k); //Jordan-Wigner string
    //         AjAN *= psidag(k);
    //         AdagjAdagN *= psi(k);
    //         AdagjAdagN *= op(sites,"F",k); //Jordan-Wigner string
    //         AdagjAdagN *= psidag(k);
    //         AdagjAN *= psi(k);
    //         AdagjAN *= op(sites,"F",k); //Jordan-Wigner string
    //         AdagjAN *= psidag(k);
    //         }

    //     AjAdagN  *= psi(N) ;
    //     AjAdagN  *= AdagN;
    //     AjAdagN  *= psidag(N);
    //     AjAN  *= psi(N) ;
    //     AjAN  *= AN;
    //     AjAN  *= psidag(N);
    //     AdagjAdagN  *= psi(N) ;
    //     AdagjAdagN  *= AdagN;
    //     AdagjAdagN  *= psidag(N);
    //     AdagjAN  *= psi(N) ;
    //     AdagjAN  *= AN;
    //     AdagjAN  *= psidag(N);
       
    //     //consider JW-transformaton the second sign is plus ,first sign is minus
    //     auto result =  elt(AjAdagN) - elt(AjAN) + elt(AdagjAN) - elt(AdagjAdagN) ;
    //     outfile10 << result << std::endl;
    //     }
    // outfile10.close();

    //calcualting  majorana <GAMMA_i gamma_2j+1>
    std::ofstream outfile9("/mnt/d/OneDrive/programs/vscodepython/marjorana/output7");
    for(int j = mid+4; j < N-1; j+=3)
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
        auto result =  -elt(A3Adagj) ;
        outfile9  << result << std::endl;
        }
    outfile9.close();

        //calcualting  majorana <GAMMA_i gamma_2j+1>
    std::ofstream outfile10("/mnt/d/OneDrive/programs/vscodepython/marjorana/output10");
    for(int j = mid+4; j < N-1; j+=3)
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
        auto result =  -elt(A3Aj) ;
        outfile10  << result << std::endl;
        }
    outfile10.close();




    //calcualting left majorana <c_j gamma_3>
    // std::ofstream outfile3("/mnt/d/OneDrive/programs/vscodepython/marjorana/output3");
    // for(int j = 6; j < N-1; j+=1)
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

    // calculate entanglment entropy for orbital A
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


    // // orbital B
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

    // // orbital C
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