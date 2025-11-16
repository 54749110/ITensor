#include "itensor/all.h"
#include "siteset_more.h"
//#include "measure.h"
using namespace itensor;
using kondoLattice = Kondo2ChannelSiteSet<ElectronSite, SpinHalfSite, ElectronSite>;
//using kondoLattice = Kondo2ChannelSiteSetv3<ElectronSite, SpinHalfSite>;
int main(int argc, char* argv[])
{
    println("//////////////////////////");
    println("Reading input file ......\n");
    //Parse the input file
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");
    //Read in individual parameters from the input file
    //second argument to getXXX methods is a default
    //in case parameter not provided in input file
    auto L = input.getInt("L");
    auto Np = input.getInt("Np");
    auto N = 3 * L;
    auto lpbc = input.getYesNo("lpbc",true);
    double t = input.getReal("t");
    double Jk = input.getReal("Jk");
    double Jh = input.getReal("Jh");
    auto readmps = input.getYesNo("readmps",false);
    auto eneropt = input.getYesNo("eneropt",true);
    auto domeas = input.getYesNo("domeas",true);
    auto quiet = input.getYesNo("quiet",true);

    // Read the sweeps parameters
    auto nsweeps = input.getInt("nsweeps");
    auto table = InputGroup(input,"sweeps");

    // suggested output file name from model parameters
    std::string runlogfile="runlog_k2c_";
    runlogfile += "L="; runlogfile += std::to_string(L);
    runlogfile += "_Np="; runlogfile += std::to_string(Np); // number of electron for each channel
    if(lpbc) runlogfile += "_PBC";
    else runlogfile += "_OBC";
    runlogfile += "_t=";  runlogfile += std::to_string(t);
    runlogfile += "_Jk="; runlogfile += std::to_string(Jk);
    runlogfile += "_Jh="; runlogfile += std::to_string(Jh);
    runlogfile += ".out";

    println("output file name suggested: ", runlogfile);

    //Create the sweeps class & print
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);


    //auto sites = kondoLattice(N);
    kondoLattice sites;
    MPS psi;
    MPO H;

    // check whether "psi_file" exists
    {
        std::ifstream ifile("psi_file");
        if(ifile) {
            readmps=true;
            println("psi_file is found, readmps is set to true"); }
        else{
            readmps=false;
            println("psi_file is NOT found, readmps is set to false"); }
    }

    if(readmps) {
        println("\n//////////////////////////////////////////////////");
        println("Reading basis, wavefunction and H from files ......");
        //kondoLattice sites;
        readFromFile("sites_file", sites);
        //MPS psi(sites);
        psi=MPS(sites);
        readFromFile("psi_file", psi);
        //MPO H(sites);
        H=MPO(sites);
        readFromFile("H_file", H);
        //auto psiHpsi = inner(psi,H,psi);
        ////auto psiHHpsi = inner(psi,H,H,psi);
        //printfln(" Intial energy information from input file: ");
        //printfln("\n<psi|H|psi> = %.10f", psiHpsi );
        ////printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
        ////printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-psiHpsi*psiHpsi );
        //println("\nTotal QN of Ground State = ",totalQN(psi));
    }
    else {
        println("\n//////////////////////////////////////////////////");
        println("Building basis, wavefunction and H from scratch ......\n");
        //
        // Initialize the site degrees of freedom.
        //
        //auto sites = kondoLattice(N);
        sites = kondoLattice(N);
        //sites = kondoLattice(N,{"ConserveNf=",true,"ConserveSz=",false});

        auto ampo = AutoMPO(sites);

        for (int j = 1; j <= (lpbc ? N : N-3); j += 3) // PBC/OBC, channel 1 hopping
        {
             ampo += -t, "Cdagup", (j+2)%N+1, "Cup", j;
             ampo += -t, "Cdagup", j,         "Cup", (j+2)%N+1;
             ampo += -t, "Cdagdn", (j+2)%N+1, "Cdn", j;
             ampo += -t, "Cdagdn", j,         "Cdn", (j+2)%N+1;
        }

        for (int j = 3; j <= (lpbc ? N : N-3); j += 3) // PBC/OBC, channel 2 hopping
        {
             ampo += -t, "Cdagup", (j+2)%N+1, "Cup", j;
             ampo += -t, "Cdagup", j,         "Cup", (j+2)%N+1;
             ampo += -t, "Cdagdn", (j+2)%N+1, "Cdn", j;
             ampo += -t, "Cdagdn", j,         "Cdn", (j+2)%N+1;
        }

        for (int j = 1; j <= N; j += 3) // channel 1 Kondo coupling
        {
            ampo += 0.5 * Jk, "S+", j, "S-", j+1;
            ampo += 0.5 * Jk, "S-", j, "S+", j+1;
            ampo +=       Jk, "Sz", j, "Sz", j+1;
        }

        for (int j = 3; j <= N; j += 3) // channel 2 Kondo coupling
        {
            ampo += 0.5 * Jk, "S+", j, "S-", j-1;
            ampo += 0.5 * Jk, "S-", j, "S+", j-1;
            ampo +=       Jk, "Sz", j, "Sz", j-1;
        }

        for (int j = 2; j <= (lpbc ? N : N-3); j += 3) // PBC/OBC, spin Heisenberg term
        {
            ampo += 0.5 * Jh, "S+", j, "S-", (j+2)%N+1;
            ampo += 0.5 * Jh, "S-", j, "S+", (j+2)%N+1;
            ampo +=       Jh, "Sz", j, "Sz", (j+2)%N+1;
        }

        H = toMPO(ampo);

        // set state
        auto state = InitState(sites);

        // channel 1
        int p = Np;
        for (int j = 1; j <= N; j += 3)
        {
            int i = N/3 - (j -1) / 3;
            if (p > i+1)
            {
                printfln("Doublely occupying site", j);

                state.set(j, "UpDn");
                p -= 2;
            }
            else if (p >0)
            {
                printfln("Singlely occupying site", j);

                state.set(j, (j%6 == 1 ? "Up" : "Dn"));
                p -= 1;
            }
            else
            {
                state.set(j, "Emp");
            }
        }

        // channel 2
        p = Np;
        for (int j = 3; j <= N; j += 3)
        {
            int i = N/3 - (j -3) / 3;
            if (p > i+1)
            {
                printfln("Doublely occupying site", j);

                state.set(j, "UpDn");
                p -= 2;
            }
            else if (p >0)
            {
                printfln("Singlely occupying site", j);

                state.set(j, (j%6 == 3 ? "Up" : "Dn"));
                p -= 1;
            }
            else
            {
                state.set(j, "Emp");
            }
        }

        // spin
        for (int j = 2; j <= N; j += 3)
        {
            if (j % 6 == 2)
                state.set(j,"Up");
            else
                state.set(j,"Dn");
        }


        // Set the initial wavefunction
        psi = MPS(state);
        //psi = randomMPS(state);

        //
        // inner calculates matrix elements of MPO's with respect to MPS's
        // inner(psi,H,psi) = <psi|H|psi>
        //
        printfln("\nInitial energy = %.5f\n", inner(psi,H,psi));
        println("\nTotal QN of Initial State = ",totalQN(psi));

    }


    if(eneropt){
        println("\n//////////////////////////////////////////////////////////////////");
        println("Beigin energy optimization, to get ground state wavefunction ......");
        //
        // Begin the DMRG calculation
        //
        auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet,"WriteDim",3100});

        // after the MPS converged, write basis, psi, and H to disk
        writeToFile("sites_file", sites);
        writeToFile("psi_file", psi);
        writeToFile("H_file", H);

        //
        // Calculate entanglement entropy
        for ( int b = 1; b <= N-1; ++b ) {
            psi.position(b);
            auto lx = leftLinkIndex(psi,b);
            auto sx = siteIndex(psi,b);
            auto [U,S,V] = svd(psi(b),{lx,sx});
            auto u = commonIndex(U,S);

            //auto wf = psi(b)*psi(b+1);
            //auto U = psi(b);
            //ITensor S, V;
            //auto spectrum = svd(wf,U,S,V);

            Real SvN = 0.;
            for ( auto n : range1(dim(u)) ) {
                auto Sn = elt(S,n,n);
                auto p = sqr(Sn);
                if ( p > 1E-12 ) SvN += -p*log(p);
            }
            printfln("Across bond b=%d, SvN = %.10f", b, SvN);
        }

        //
        // Print the final energy reported by DMRG
        //
        println("\nTotal QN of Ground State = ",totalQN(psi));
        printfln("\nGround State Energy = %.10f", energy);
        printfln("\n<psi|H|psi> / N = %.10f", energy/N );

    }

    if(domeas){
        auto ampo_ek = AutoMPO(sites);
        for (int j = 1; j <= (lpbc ? N : N-3); j += 3) // PBC/OBC, channel 1 hopping
        {
             ampo_ek += -t, "Cdagup", (j+2)%N+1, "Cup", j;
             ampo_ek += -t, "Cdagup", j,         "Cup", (j+2)%N+1;
             ampo_ek += -t, "Cdagdn", (j+2)%N+1, "Cdn", j;
             ampo_ek += -t, "Cdagdn", j,         "Cdn", (j+2)%N+1;
        }

        for (int j = 3; j <= (lpbc ? N : N-3); j += 3) // PBC/OBC, channel 2 hopping
        {
             ampo_ek += -t, "Cdagup", (j+2)%N+1, "Cup", j;
             ampo_ek += -t, "Cdagup", j,         "Cup", (j+2)%N+1;
             ampo_ek += -t, "Cdagdn", (j+2)%N+1, "Cdn", j;
             ampo_ek += -t, "Cdagdn", j,         "Cdn", (j+2)%N+1;
        }

        auto ampo_ejk = AutoMPO(sites);
        for (int j = 1; j <= N; j += 3) // channel 1 Kondo coupling
        {
            ampo_ejk += 0.5 * Jk, "S+", j, "S-", j+1;
            ampo_ejk += 0.5 * Jk, "S-", j, "S+", j+1;
            ampo_ejk +=       Jk, "Sz", j, "Sz", j+1;
        }

        for (int j = 3; j <= N; j += 3) // channel 2 Kondo coupling
        {
            ampo_ejk += 0.5 * Jk, "S+", j, "S-", j-1;
            ampo_ejk += 0.5 * Jk, "S-", j, "S+", j-1;
            ampo_ejk +=       Jk, "Sz", j, "Sz", j-1;
        }

        auto ampo_ejh = AutoMPO(sites);
        for (int j = 2; j <= (lpbc ? N : N-3); j += 3) // PBC/OBC, spin Heisenberg term
        {
            ampo_ejh += 0.5 * Jh, "S+", j, "S-", (j+2)%N+1;
            ampo_ejh += 0.5 * Jh, "S-", j, "S+", (j+2)%N+1;
            ampo_ejh +=       Jh, "Sz", j, "Sz", (j+2)%N+1;
        }

        auto hk = toMPO(ampo_ek);
        auto hjk = toMPO(ampo_ejk);
        auto hjh = toMPO(ampo_ejh);
        printfln("\nKinetic Energy = %.10f", inner(psi,hk,psi));
        printfln("\nKondo part Energy = %.10f", inner(psi,hjk,psi));
        printfln("\nHeisenberg part Energy = %.10f", inner(psi,hjh,psi));

        println("\n////////////////////////////");
        println("Start to perform measurement of spin correlation\n");
        //
        // Measure Si.Sj of every {i,j}, and total M
        //
        auto totalM = 0.0;
        auto Msquare = 0.0;
        auto Mzsquare = 0.0;
        std::vector<double> SiSj_meas={};
        std::vector<double> SiSjzz_meas={};
        std::vector<double> SiSjpm_meas={};
        std::vector<double> Sz_meas={};
        std::vector<Cplx> Sp_meas={};
        std::vector<Cplx> Sm_meas={};
        for ( int i = 2; i <= N; i += 3 ) {
            //'gauge' the MPS to site i
            psi.position(i); 
            
            //psinc(1) *= psi(0); //Uncomment if doing iDMRG calculation

            // i == j part

            // magnetization
            auto ket = psi(i);
            auto bra = dag(prime(ket,"Site"));
            auto sz_tmp = ((bra*sites.op("Sz",i)*ket).cplx()).real();
            Sz_meas.emplace_back(sz_tmp);
            totalM +=  sz_tmp;

            auto sp_tmp = (bra*sites.op("S+",i)*ket).cplx();
            Sp_meas.emplace_back(sp_tmp);
            auto sm_tmp = (bra*sites.op("S-",i)*ket).cplx();
            Sm_meas.emplace_back(sm_tmp);

            auto szsz_tmp = 0.0;
            szsz_tmp = (( prime(bra*sites.op("Sz",i),"Site")*sites.op("Sz",i)*ket).cplx()).real();
            SiSjzz_meas.emplace_back(szsz_tmp);
            Mzsquare += szsz_tmp;
            auto spsm_tmp = 0.0;
            spsm_tmp = (( prime(bra*sites.op("S+",i),"Site")*sites.op("S-",i)*ket).cplx()).real();
            SiSjpm_meas.emplace_back(spsm_tmp);

            auto ss_tmp = szsz_tmp + spsm_tmp;
            Msquare += ss_tmp;
            SiSj_meas.emplace_back(ss_tmp);
            println( i, " ", i, " ", ss_tmp );
            
            if ( i < N ) {
                // i != j part
                //index linking i to i+1:
                auto ir = commonIndex(psi(i),psi(i+1),"Link");
   
                auto op_ip = sites.op("S+",i);
                auto op_im = sites.op("S-",i);
                auto op_iz = sites.op("Sz",i);
                auto Cpm = psi(i)*op_ip*dag(prime(prime(psi(i),"Site"),ir));
                auto Cmp = psi(i)*op_im*dag(prime(prime(psi(i),"Site"),ir));
                auto Czz = psi(i)*op_iz*dag(prime(prime(psi(i),"Site"),ir));
                for(int j = i+1; j <= N; ++j) {
                    Cpm *= psi(j);
                    Cmp *= psi(j);
                    Czz *= psi(j);

                    if( j%3 == 2 ){

                        auto jl = commonIndex(psi(j),psi(j-1),"Link");

                        auto op_jm = sites.op("S-",j);
                        auto ss_tmp = 0.0;
                        ss_tmp += 0.5*(( (Cpm*op_jm)*dag(prime(prime(psi(j),jl),"Site")) ).cplx()).real();

                        auto op_jp = sites.op("S+",j);
                        ss_tmp += 0.5*(( (Cmp*op_jp)*dag(prime(prime(psi(j),jl),"Site")) ).cplx()).real();

                        auto spm_tmp = ss_tmp;
                        SiSjpm_meas.emplace_back(ss_tmp);

                        auto op_jz = sites.op("Sz",j);
                        ss_tmp += (( (Czz*op_jz)*dag(prime(prime(psi(j),jl),"Site")) ).cplx()).real();

                        SiSjzz_meas.emplace_back(ss_tmp-spm_tmp);
                        Mzsquare += (ss_tmp - spm_tmp)*2.0;
                        
                        SiSj_meas.emplace_back(ss_tmp);
                        Msquare += ss_tmp*2.0;
                        println( i, " ", j, " ", ss_tmp ); 
                    }

                    if(j < N) {
                        Cpm *= dag(prime(psi(j),"Link"));
                        Cmp *= dag(prime(psi(j),"Link"));
                        Czz *= dag(prime(psi(j),"Link"));
                    }
                }
            }
        }
        printfln("Total M = %.10e", totalM );
        printfln("Msquare = %.10e", Msquare );
        printfln("Mzsquare = %.10e", Mzsquare );

        std::ofstream fSzout("Siz.out",std::ios::out);
        fSzout.precision(16);
        for (std::vector<double>::const_iterator i = Sz_meas.begin(); i != Sz_meas.end(); ++i)
                fSzout << *i << ' ';

        std::ofstream fSpout("Sip.out",std::ios::out);
        fSpout.precision(16);
        for (std::vector<Cplx>::const_iterator i = Sp_meas.begin(); i != Sp_meas.end(); ++i)
                fSpout << *i << ' ';

        std::ofstream fSmout("Sim.out",std::ios::out);
        fSmout.precision(16);
        for (std::vector<Cplx>::const_iterator i = Sm_meas.begin(); i != Sm_meas.end(); ++i)
                fSmout << *i << ' ';

        std::ofstream fSiSjout("SiSj.out",std::ios::out);
        fSiSjout.precision(16);
        for (std::vector<double>::const_iterator i = SiSj_meas.begin(); i != SiSj_meas.end(); ++i)
                fSiSjout << *i << ' ';

        std::ofstream fSiSjzzout("SiSjzz.out",std::ios::out);
        fSiSjzzout.precision(16);
        for (std::vector<double>::const_iterator i = SiSjzz_meas.begin(); i != SiSjzz_meas.end(); ++i)
                fSiSjzzout << *i << ' ';

        std::ofstream fSiSjpmout("SiSjpm.out",std::ios::out);
        fSiSjpmout.precision(16);
        for (std::vector<double>::const_iterator i = SiSjpm_meas.begin(); i != SiSjpm_meas.end(); ++i)
                fSiSjpmout << *i << ' ';

        println("\n////////////////////////////");
        println("Start to perform measurement of fermion two body correlation\n");
        //
        // Measure cj^dag.ci of every {i,j}, for j>=i
        //
        std::vector<double> cicjpm1_meas={};
        std::vector<double> cicjpm2_meas={};
        // Note we are applying operator to right (to on the ket)
        // Also note the basis has larger index on the right side,
        // and we apply operator with larger index first to right (to on the ket)
        // calculate <c_j^dag c_i> for j>i
        // calculate <c_j^dag c_i> = -(a_iup F_i) F_i+1 F_i+2 ... F_j-1 a_jup^dag - a_idn F_i+1 F_i+2 ... F_j-1 (F_j a_jdn^dag)
        for ( int j = N-2; j >= 1; j -= 3 ) {
            //'gauge' the MPS to site j
            psi.position(j); 
            
            //psinc(1) *= psi(0); //Uncomment if doing iDMRG calculation

            // j == i part

            auto ket = psi(j);
            auto bra = dag(prime(ket,"Site"));

            auto cpcm_tmp = ((bra*sites.op("Nup",j)*ket).cplx()).real() + ((bra*sites.op("Ndn",j)*ket).cplx()).real();
            cicjpm1_meas.emplace_back(cpcm_tmp);
            println( j, " ", j, " ", cpcm_tmp );

            if ( j > 1 ) {
                // i != j part
                //index linking i to i+1:
                auto jr = commonIndex(psi(j),psi(j-1),"Link");
   
                auto op_jup = sites.op("Adagup",j);
                auto op_jdn = sites.op("Adagdn",j);
                // Note we are applying operator to right (to on the ket)
                // Also note the basis has larger index on the right side,
                // and we apply operator with larger index first to right (to on the ket)
                auto Cpmup = psi(j)*op_jup*dag(prime(prime(psi(j),"Site"),jr));
                auto Cpmdn = noPrime(psi(j)*op_jdn,"Site")*sites.op("F",j)*dag(prime(prime(psi(j),"Site"),jr));

                for(int i = j-1; i >= 1; --i) {

                    if( i%3 == 0 ) {
                        // fermion site
                        Cpmup = noPrime(Cpmup*psi(i)*sites.op("F",i),"Site");
                        Cpmdn = noPrime(Cpmdn*psi(i)*sites.op("F",i),"Site");
                    }
                    else if( i%3 == 2) {
                        Cpmup *= psi(i);
                        Cpmdn *= psi(i);
                    }
                    else {
                        Cpmup *= psi(i);
                        Cpmdn *= psi(i);

                        auto il = commonIndex(psi(i),psi(i+1),"Link");
                        auto op_iup = sites.op("Aup",i);
                        auto op_idn = sites.op("Adn",i);
                        auto cc_tmp = 0.0;
                        cc_tmp += -(( noPrime(Cpmup*sites.op("F",i),"Site")*op_iup*dag(prime(prime(psi(i),il),"Site")) ).cplx()).real();
                        cc_tmp += -(( (Cpmdn*op_idn)*dag(prime(prime(psi(i),il),"Site")) ).cplx()).real();
                        cicjpm1_meas.emplace_back(cc_tmp);

                        println( i, " ", j, " ", cc_tmp ); 

                        Cpmup = noPrime(Cpmup*sites.op("F",i),"Site");
                        Cpmdn = noPrime(Cpmdn*sites.op("F",i),"Site");
                    }

                    Cpmup *= dag(prime(psi(i),"Link"));
                    Cpmdn *= dag(prime(psi(i),"Link"));
                }
            }
        }

        for ( int j = N; j >= 1; j -= 3 ) {
            //'gauge' the MPS to site j
            psi.position(j); 
            
            //psinc(1) *= psi(0); //Uncomment if doing iDMRG calculation

            // j == i part

            auto ket = psi(j);
            auto bra = dag(prime(ket,"Site"));

            auto cpcm_tmp = ((bra*sites.op("Nup",j)*ket).cplx()).real() + ((bra*sites.op("Ndn",j)*ket).cplx()).real();
            cicjpm2_meas.emplace_back(cpcm_tmp);
            println( j, " ", j, " ", cpcm_tmp );

            if ( j > 1 ) {
                // i != j part
                //index linking i to i+1:
                auto jr = commonIndex(psi(j),psi(j-1),"Link");
   
                auto op_jup = sites.op("Adagup",j);
                auto op_jdn = sites.op("Adagdn",j);
                // Note we are applying operator to right (to on the ket)
                // Also note the basis has larger index on the right side,
                // and we apply operator with larger index first to right (to on the ket)
                auto Cpmup = psi(j)*op_jup*dag(prime(prime(psi(j),"Site"),jr));
                auto Cpmdn = noPrime(psi(j)*op_jdn,"Site")*sites.op("F",j)*dag(prime(prime(psi(j),"Site"),jr));

                for(int i = j-1; i >= 1; --i) {

                    if( i%3 == 1 ) {
                        // fermion site
                        Cpmup = noPrime(Cpmup*psi(i)*sites.op("F",i),"Site");
                        Cpmdn = noPrime(Cpmdn*psi(i)*sites.op("F",i),"Site");
                    }
                    else if( i%3 == 2) {
                        Cpmup *= psi(i);
                        Cpmdn *= psi(i);
                    }
                    else {
                        Cpmup *= psi(i);
                        Cpmdn *= psi(i);

                        auto il = commonIndex(psi(i),psi(i+1),"Link");
                        auto op_iup = sites.op("Aup",i);
                        auto op_idn = sites.op("Adn",i);
                        auto cc_tmp = 0.0;
                        cc_tmp += -(( noPrime(Cpmup*sites.op("F",i),"Site")*op_iup*dag(prime(prime(psi(i),il),"Site")) ).cplx()).real();
                        cc_tmp += -(( (Cpmdn*op_idn)*dag(prime(prime(psi(i),il),"Site")) ).cplx()).real();
                        cicjpm2_meas.emplace_back(cc_tmp);

                        println( i, " ", j, " ", cc_tmp ); 

                        Cpmup = noPrime(Cpmup*sites.op("F",i),"Site");
                        Cpmdn = noPrime(Cpmdn*sites.op("F",i),"Site");
                    }

                    Cpmup *= dag(prime(psi(i),"Link"));
                    Cpmdn *= dag(prime(psi(i),"Link"));
                }
            }
        }

        std::ofstream fcicjpm1out("cicjpm1.out",std::ios::out);
        fcicjpm1out.precision(16);
        for (std::vector<double>::const_iterator i = cicjpm1_meas.begin(); i != cicjpm1_meas.end(); ++i)
                fcicjpm1out << *i << ' ';

        std::ofstream fcicjpm2out("cicjpm2.out",std::ios::out);
        fcicjpm2out.precision(16);
        for (std::vector<double>::const_iterator i = cicjpm2_meas.begin(); i != cicjpm2_meas.end(); ++i)
                fcicjpm2out << *i << ' ';


        println("\n////////////////////////////");
        println("Start to perform measurement of bb correlation\n");
        //
        // Measure cdni.S^-(i).S^+(i+1)....S^-(j).cupj^dag of every {i,j}, for j>=i, if j-i+1 is odd
        //      or cupi.S^+(i).S^-(i+1)....S^-(j).cupj^dag of every {i,j}, for j>=i, if j-i+1 is even
        //
        // <cdni cupj^dag> = - a_idn      F_i+1 F_i+2 ... F_j-1 a_jup^dag
        // <cupi cupj^dag> =  (a_iup F_i) F_i+1 F_i+2 ... F_j-1 a_jup^dag
        std::vector<double> bibj1_meas={};
        // Note we are applying operator to right (to on the ket)
        // Also note the basis has larger index on the right side,
        // and we apply operator with larger index first to right (to on the ket)
        for ( int j = N-1; j >= 1; j -= 3 ) {
            //'gauge' the MPS to site j
            psi.position(j); 
            
            //psinc(1) *= psi(0); //Uncomment if doing iDMRG calculation

            // j == i part

            if ( j > 1 ) {
                //index linking j to j-1:
                auto jr = commonIndex(psi(j),psi(j-1),"Link");
   
                // Note we are applying operator to right (to on the ket)
                // Also note the basis has larger index on the right side,
                // and we apply operator with larger index first to right (to on the ket)
                // current site spin is on the right hand side of channel-1 fermion
                auto bibjop = psi(j)*sites.op("S-",j)*dag(prime(prime(psi(j),"Site"),jr));
                bibjop = noPrime(bibjop*psi(j-1)*sites.op("Adagup",j-1),"Site");
                auto bb_tmp = -((bibjop*sites.op("Adn",j-1)*dag(prime(prime(psi(j-1),"Site"),jr))).cplx()).real();
                bibj1_meas.emplace_back(bb_tmp);
                println( j-1, " ", j-1, " ", bb_tmp );

                bibjop *= dag(prime(psi(j-1),"Link"));

                for(int i = j-2; i >= 1; --i) {

                    if( i%3 == 0 ) {
                        // fermion site
                        bibjop = noPrime(bibjop*psi(i)*sites.op("F",i),"Site");
                    }
                    else if( i%3 == 2) {
                        // spin site
                        if( (j-i+1)%6 == 1 ) {
                            bibjop = noPrime(bibjop*psi(i)*sites.op("S-",i),"Site");
                        }
                        else {
                            bibjop = noPrime(bibjop*psi(i)*sites.op("S+",i),"Site");
                        }
                    }
                    else {
                        bibjop *= psi(i);

                        auto il = commonIndex(psi(i),psi(i+1),"Link");
                        auto bb_tmp = 0.0;
                        if( (j-i+1)%6 == 2 ) {
                            auto op_idn = sites.op("Adn",i);
                            bb_tmp = -(( (bibjop*op_idn)*dag(prime(prime(psi(i),il),"Site")) ).cplx()).real();
                        }
                        else {
                            auto op_iup = sites.op("Aup",i);
                            bb_tmp = (( noPrime(bibjop*sites.op("F",i),"Site")*op_iup*dag(prime(prime(psi(i),il),"Site")) ).cplx()).real();
                        }

                        bibj1_meas.emplace_back(bb_tmp);
                        println( i, " ", j-1, " ", bb_tmp ); 

                        bibjop = noPrime(bibjop*sites.op("F",i),"Site");
                    }

                    bibjop *= dag(prime(psi(i),"Link"));
                }
            }
        }

        std::vector<double> bibj2_meas={};
        for ( int j = N; j >= 1; j -= 3 ) {
            //'gauge' the MPS to site j
            psi.position(j); 
            
            //psinc(1) *= psi(0); //Uncomment if doing iDMRG calculation

            // j == i part

            if ( j > 1 ) {
                //index linking j to j-1:
                auto jr = commonIndex(psi(j),psi(j-1),"Link");
   
                // Note we are applying operator to right (to on the ket)
                // Also note the basis has larger index on the right side,
                // and we apply operator with larger index first to right (to on the ket)
                auto bbop_tmp = noPrime(psi(j)*sites.op("Adagup",j),"Site")*sites.op("Adn",j)*dag(prime(prime(psi(j),"Site"),jr));
                // current site spin is on the left hand side of channel-2 fermion
                auto bb_tmp = -((noPrime(bbop_tmp*psi(j-1)*sites.op("S-",j-1),"Site")*dag(prime(psi(j-1),jr))).cplx()).real();
                bibj2_meas.emplace_back(bb_tmp);
                println( j, " ", j, " ", bb_tmp );

                auto bibjop = psi(j)*sites.op("Adagup",j)*dag(prime(prime(psi(j),"Site"),jr));
                bibjop = noPrime(bibjop*psi(j-1)*sites.op("S-",j-1),"Site");
                bibjop *= dag(prime(psi(j-1),"Link"));

                for(int i = j-2; i >= 1; --i) {

                    if( i%3 == 1 ) {
                        // fermion site
                        bibjop = noPrime(bibjop*psi(i)*sites.op("F",i),"Site");
                    }
                    else if( i%3 == 2) {
                        // spin site
                        if( (j-i+1)%6 == 2 ) {
                            bibjop = noPrime(bibjop*psi(i)*sites.op("S-",i),"Site");
                        }
                        else {
                            bibjop = noPrime(bibjop*psi(i)*sites.op("S+",i),"Site");
                        }
                    }
                    else {
                        bibjop *= psi(i);

                        auto bb_tmp = 0.0;
                        if( (j-i+1)%6 == 1 ) {
                            auto op_idn = sites.op("Adn",i);
                            auto bbop_tmp = bibjop*op_idn*dag(prime(prime(psi(i),"Link"),"Site"));
                            auto il = commonIndex(psi(i),psi(i-1),"Link");
                            // current site spin is on the left hand side of channel-2 fermion
                            bb_tmp = -((bbop_tmp*psi(i-1)*sites.op("S-",i-1)*dag(prime(prime(psi(i-1),il),"Site")) ).cplx()).real();
                        }
                        else {
                            auto op_iup = sites.op("Aup",i);
                            auto bbop_tmp = noPrime(bibjop*sites.op("F",i),"Site")*op_iup*dag(prime(prime(psi(i),"Link"),"Site"));
                            auto il = commonIndex(psi(i),psi(i-1),"Link");
                            // current site spin is on the left hand side of channel-2 fermion
                            bb_tmp = ((bbop_tmp*psi(i-1)*sites.op("S+",i-1)*dag(prime(prime(psi(i-1),il),"Site")) ).cplx()).real();
                        }

                        bibj2_meas.emplace_back(bb_tmp);
                        println( i, " ", j, " ", bb_tmp ); 

                        bibjop = noPrime(bibjop*sites.op("F",i),"Site");
                    }

                    bibjop *= dag(prime(psi(i),"Link"));
                }
            }
        }
        std::ofstream fbibj1out("bibj1.out",std::ios::out);
        fbibj1out.precision(16);
        for (std::vector<double>::const_iterator i = bibj1_meas.begin(); i != bibj1_meas.end(); ++i)
                fbibj1out << *i << ' ';

        std::ofstream fbibj2out("bibj2.out",std::ios::out);
        fbibj2out.precision(16);
        for (std::vector<double>::const_iterator i = bibj2_meas.begin(); i != bibj2_meas.end(); ++i)
                fbibj2out << *i << ' ';
    }

    return 0;
}
