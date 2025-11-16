#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int 
main()
    {
    // N must be 3Z
    int N = 60;
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

    // kitaev model 在 lieb lattice 上面 ，由于 kitaev model p-wave 超导，自旋同向。
    // 则这里默认自旋向上，每个格点只能空或是自旋上的一个电子。可以选取FermionSite。
    auto sites = Fermion(N,{"ConserveQNs=", false}); 


    auto t= 0.0003;
    auto J= 1;
    auto V= 0.2;
    auto mu = 0; 
    auto D =0.6*t;
    // auto t=1;
    // auto J=1;
    // auto V=1;
    // auto mu = 0; 
    // auto D =1 ;

    // Open boundary condition
    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-3 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+3, "C", j;
            ampo += -t, "Cdag", j,   "C", j+3;
    }
    for(int j = 1; j <= N-5; j += 3) // electron hopping over orbitals
        {
            ampo += -J, "Cdag", j+1, "C", j;
            ampo += -J, "Cdag", j,   "C", j+1;
            ampo += -J, "Cdag", j+2, "C", j;
            ampo += -J, "Cdag", j,   "C", j+2;
            ampo += -J, "Cdag", j+4, "C", j;
            ampo += -J, "Cdag", j,   "C", j+4;
            ampo += -J, "Cdag", j+5, "C", j;
            ampo += -J, "Cdag", j,   "C", j+5;
    }
        
    for (int j = 2; j <= N-1 ; j += 3) //onsite V 这里几项我尝试配做particle-hole对称的形式，去掉了一个常数项
        {
            ampo += V, "Cdag", j, "C", j;
          //  ampo += -0.5*V, "C", j, "Cdag", j;
    }

    for (int j = 3; j <= N ; j += 3) //onsite V
        {
            ampo += -V, "Cdag", j, "C", j;
          //  ampo +=  0.5*V, "C", j, "Cdag", j;
    }

    for (int j = 1; j <= N ; j += 1) //chemical potential
        {
            // ampo += -0.5*mu-0.0001, "Cdag", j, "C", j;
            // ampo += 0.5*mu-0.0001, "C", j, "Cdag", j;
            ampo += -0.5*mu, "Cdag", j, "C", j;
            ampo += 0.5*mu, "C", j, "Cdag", j;
    }

    for (int j = 1; j <= N-3 ; j += 1) //electron pairing
        {
            ampo += D, "Cdag", j+3, "Cdag", j;
            ampo += D, "C", j,   "C", j+3;
    }
    
    ampo += -J, "Cdag", N-1, "C", N-2;
    ampo += -J, "Cdag", N-2,   "C", N-1;
    ampo += -J, "Cdag", N, "C", N-2;
    ampo += -J, "Cdag", N-2,   "C", N;

    H = toMPO(ampo);
    // Several tries to the occupation of initial state
    //
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%3 == 0) state.set(i,"Occ");
        else         state.set(i,"Emp");
        }
    psi0 = MPS(state);
    
    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //

     

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(80);
    sweeps.maxdim() = 10,10,10,10,20,20,20,20,40,40,40,40,100,100,100,100,200,200,200,200,400,400,400,400,600,600,600,600;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = -1E-3,-1E-4,-1E-5,-1E-6,-1E-7,-1E-9,0.0;
    //sweeps.noise() = 1E-1,1E-1,1E-2,1E-2,1E-2,1E-2,1E-2,1E-2,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3,1E-4,1E-4,1E-4,1E-4,1E-4,1E-4,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-9,1E-9,1E-9,1E-9,1E-9,1E-9,1E-10,1E-10,1E-10,1E-10,1E-10,1E-10,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,0.0;
    //println(sweeps);

    //
    // Begin the DMRG calculation
    //
   auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
   printfln("Initial energy = %.5f", inner(psi0,H,psi0) );
   printfln("\nGround State Energy = %.10f",energy);

   for(int j = 1; j <= N; ++j)
        {
        //re-gauge psi to get ready to measure at position j
        auto ampo3 = AutoMPO(sites) ;
        ampo3 += 1, "Cdag", j, "C", 1;
        G1 = toMPO(ampo3);
        auto ampo4 = AutoMPO(sites) ;
        ampo4 += 1, "Cdag", j, "Cdag", 1;
        G2 = toMPO(ampo4);
        auto ampo5 = AutoMPO(sites) ;
        ampo5 += -1, "Cdag", j, "C", N;
        G3 = toMPO(ampo5);
        auto ampo6 = AutoMPO(sites) ;
        ampo6 += 1, "Cdag", j, "Cdag", N;
        G4 = toMPO(ampo6);
        auto wave1 = inner(psi , G1 , psi) ;
        auto wave2 = inner(psi , G2 , psi) ;
        auto wave3 = inner(psi , G3 , psi) ;
        auto wave4 = inner(psi , G4 , psi) ;
        auto majorana = wave1 + wave2 + wave3 + wave4 ;
        printfln("Site %d occupation: %.5f", j, majorana);
    }

    // PBC
    auto ampo1 = AutoMPO(sites);
    for (int j = 1; j <= N-3 ; j += 1) //electron hopping 
        {
            ampo1 += -t, "Cdag", j+3, "C", j;
            ampo1 += -t, "Cdag", j,   "C", j+3;
    }
    for(int j = 1; j <= N-5; j += 3) // electron hopping over orbitals
        {
            ampo1 += -J, "Cdag", j+1, "C", j;
            ampo1 += -J, "Cdag", j,   "C", j+1;
            ampo1 += -J, "Cdag", j+2, "C", j;
            ampo1 += -J, "Cdag", j,   "C", j+2;
            ampo1 += -J, "Cdag", j+4, "C", j;
            ampo1 += -J, "Cdag", j,   "C", j+4;
            ampo1 += -J, "Cdag", j+5, "C", j;
            ampo1 += -J, "Cdag", j,   "C", j+5;
    }
        
    for (int j = 2; j <= N-1 ; j += 3) //onsite V 这里几项我尝试配做particle-hole对称的形式，去掉了一个常数项
        {
            ampo1 += 0.5*V, "Cdag", j, "C", j;
            ampo1 += -0.5*V, "C", j, "Cdag", j;
    }

    for (int j = 3; j <= N ; j += 3) //onsite V
        {
            ampo1 += -0.5*V, "Cdag", j, "C", j;
            ampo1 +=  0.5*V, "C", j, "Cdag", j;
    }

    for (int j = 1; j <= N ; j += 1) //chemical potential
        {
            ampo1 += -0.5*mu-0.0001, "Cdag", j, "C", j;
            ampo1 += 0.5*mu-0.0001, "C", j, "Cdag", j;
            // ampo += -0.5*mu, "Cdag", j, "C", j;
            // ampo += 0.5*mu, "C", j, "Cdag", j;
    }

    for (int j = 1; j <= N-3 ; j += 1) //electron pairing
        {
            ampo1 += D, "Cdag", j+3, "Cdag", j;
            ampo1 += D, "C", j,   "C", j+3;
    }
    
    ampo1 += -J, "Cdag", N-1, "C", N-2;
    ampo1 += -J, "Cdag", N-2,   "C", N-1;
    ampo1 += -J, "Cdag", N, "C", N-2;
    ampo1 += -J, "Cdag", N-2,   "C", N;
    

    // PBC condition
    ampo1 += -t, "Cdag", N, "C", 3;
    ampo1 += -t, "Cdag", 3,   "C", N;
    ampo1 += -t, "Cdag", N-1, "C", 2;
    ampo1 += -t, "Cdag", 2,   "C", N-1;
    ampo1 += -t, "Cdag", N-2, "C", 1;
    ampo1 += -t, "Cdag", 1,   "C", N-2;
    ampo1 += D, "Cdag", N, "Cdag", 3;
    ampo1 += D, "C", 3,   "C", N;
    ampo1 += D, "Cdag", N-1, "Cdag", 2;
    ampo1 += D, "C", 2,   "C", N-1;
    ampo1 += D, "Cdag", N-2, "Cdag", 1;
    ampo1 += D, "C", 1,   "C", N-2;
    ampo1 += -J, "Cdag", 2, "C", N-2;
    ampo1 += -J, "Cdag", N-2,   "C", 2;
    ampo1 += -J, "Cdag", 3, "C", N-2;
    ampo1 += -J, "Cdag", N-2,   "C", 3;


    Hpbc = toMPO(ampo1);

    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    //printfln("Initial energy = %.5f", inner(psi0,Hpbc,psi0) );
    
    auto sweeps2 = Sweeps(20);
    sweeps2.maxdim() = 2,2,2,2,4,4,4,4,8,8,8,8,16,16,16,16,32,32,32,32; //640,1280 ;
    sweeps2.cutoff() = 1E-9;
    sweeps2.niter() = 2;
    sweeps2.noise() = 1E-7,1E-8,1E-9,1E-10,0.0;
    //
    // Begin the DMRG calculation
    //
   //auto [energypbc,psipbc] = dmrg(Hpbc,psi0,sweeps2,"Quiet");
  // printfln("\nGround State pbc Energy = %.10f",energypbc);


    
    // 利用投影的方法作用到psi_obc得到majorana 波函数， 投影算符为 P = I-|psi_pbc><psi_pbc|
    // for(int j = 1; j <= N; ++j)
    //     {
    //     //re-gauge psi to get ready to measure at position j
    //     auto ampo2 = AutoMPO(sites) ;
    //     ampo2 += 1, "Cdag", j, "C", j;
    //     n_j = toMPO(ampo2);
    //     auto wave1 = inner(psi , n_j , psi) ;
    //     auto wave2 = inner(psipbc , n_j , psi) ;
    //     auto wave3 = inner(psi , n_j , psipbc) ;
    //     auto wave4 = inner(psipbc , n_j , psipbc) ;
    //     auto b     = inner(psipbc,psi) ;
    //     auto majorana = wave1 - b*wave2 - b*wave3 + b*b*wave4 ;
    //     printfln("Site %d occupation: %.5f", j, majorana);
    // }
    


    //计算第一激发态，基态再填充一个electron 正好填到majorana能级
   // auto wfs = std::vector<MPS>(1);
   // wfs.at(0) = psi;

    // Here the Weight option sets the energy penalty for
    // psi1 having any overlap with psi
    //auto [en1,psi1] = dmrg(H,wfs,psi0,sweeps,"Quiet");

    //printfln("\nExcited State Energy = %.10f",en1);
    
    // 马约拉纳波函数 =<基态|c^{+}_i|第一激发态>+<第一激发态|c^{+}_i|基态>
    // 参考文献：https://arxiv.org/pdf/1104.5493
    // for(int j = 1; j <= N; ++j)
    //    {
    //    psi.position(j);
    //    psi1.position(j);
    //    auto ket0 = psi(j);
    //    auto bra0 = dag(prime(ket0,"Site")).noPrime("Link");;
    //    auto ket1 = psi1(j);
    //    auto bra1 = dag(prime(ket1,"Site")).noPrime("Link");;
    //    auto c_i= op(sites, "Adag", j) ;
    //    auto wave1 = elt(bra1*c_i*ket0) ;
    //    auto wave2 = elt(bra0*c_i*ket1) ;
    //    auto majorana = wave1 + wave2 ;
    //    printfln("Site %d occupation: %.5f", j, majorana);
    // }


    // // itensor 输出的psi 是一个MPS 但是可以通过计算 n_i 的期望值得到近似特征向量。
    // for(int j = 1; j <= N; ++j)
    //     {
    //     //re-gauge psi to get ready to measure at position j
    //     psi.position(j);
    //     auto ket = psi(j);
    //     auto bra = dag(prime(ket,"Site"));
    //     auto n_j = op(sites, "N", j); // 平均占据数
    //     auto wave = elt(bra*n_j*ket); 
    //     printfln("Site %d occupation: %.5f", j, wave);
    // }
    

    // for(int j = 1; j <= N; ++j)
    //     {
    //    auto ampo1 = AutoMPO(sites);
    //    ampo1 += 1, "Adag", j ;

    //    H1 = toMPO(ampo1);
    //    auto wave1 = inner(psi1,H1,psi) ;
    //    auto wave2 = inner(psi,H1,psi1) ;
    //    auto majorana1 = wave1 + wave2 ;
    //    auto majorana2 = wave1 - wave2 ;
    //    printfln("Site %d occupation 1: %.5f", j, majorana1);
    //    printfln("Site %d occupation 2: %.5f", j, majorana2);
    // }

    return 0;
    }