#include "itensor/all.h"
#include "itensor/util/print_macro.h"
using namespace itensor;


// std::vector<double> QML(const std::vector<double>& inputArray) {
//     if (inputArray.size() < 2) {
//         return {}; // 如果输入数组太小，返回空数组
//     }

//     std::vector<double> absLogArray;
//     for (double num : inputArray) {
//         absLogArray.push_back(log(abs(num))); // 先取绝对值，再取对数
//     }

//     std::vector<double> result;
//     for (size_t i = 1; i < absLogArray.size(); ++i) {
//         result.push_back(absLogArray[i] - absLogArray[i-1]); // 计算相邻差值
//     }

//     return result; // 返回差值数组
// }


int 
main()
    {
    int N = 40;
    MPO H;
    MPS psi0;
    MPO H1 ;
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

    auto miu = 0.1; 
    auto t=1;
    auto D=1;

    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N-1 ; j += 1) //electron hopping 
        {
            ampo += -t, "Cdag", j+1, "C", j;
            ampo += -t, "Cdag", j,   "C", j+1;
        
    }
    for (int j = 1; j <= N-1 ; j += 1) //electron pairing
        {
           ampo += D, "Cdag", j+1, "Cdag", j;
           ampo += D, "C", j,   "C", j+1;
    }

    for (int j = 1; j <= N ; j += 1) //chemical potential
    {
        // ampo += -0.5*mu-0.0001, "Cdag", j, "C", j;
        // ampo += 0.5*mu-0.0001, "C", j, "Cdag", j;
        ampo += -0.5*miu, "Cdag", j, "C", j;
        ampo += 0.5*miu, "C", j, "Cdag", j;
    }

    // auto ampopbc = AutoMPO(sites);
    // for (int j = 1; j <= N-1 ; j += 1) //electron hopping 
    //     {
    //         ampopbc += -t, "Cdag", j+1, "C", j;
    //         ampopbc += -t, "Cdag", j,   "C", j+1;
        
    // }
    // for (int j = 1; j <= N-1 ; j += 1) //electron pairing
    //     {
    //        ampopbc += D, "Cdag", j+1, "Cdag", j;
    //        ampopbc += D, "C", j,   "C", j+1;
    // }

    // for (int j = 1; j <= N ; j += 1) //chemical potential
    // {
    //     // ampo += -0.5*mu-0.0001, "Cdag", j, "C", j;
    //     // ampo += 0.5*mu-0.0001, "C", j, "Cdag", j;
    //     ampopbc += -0.5*miu, "Cdag", j, "C", j;
    //     ampopbc += 0.5*miu, "C", j, "Cdag", j;
    // }
    // // ampopbc += D, "Cdag", 1, "Cdag", N;
    // // ampopbc+= D, "C", N,   "C", 1;
    // ampopbc += D, "Cdag", N, "Cdag", 1;
    // ampopbc+= D, "C", 1,   "C", N;
    // ampopbc += -t, "Cdag", N, "C", 1;
    // ampopbc += -t, "Cdag", 1,   "C", N;

    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%2 == 1) state.set(i,"Emp");
        else         state.set(i,"Occ");
        }
    psi0 = MPS(state);


    H = toMPO(ampo);
    // H1 = toMPO(ampopbc);
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
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10,20,100,100,200 ;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,0.0;
    //println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
    //auto [energypbc,psi1] = dmrg(H1,psi0,sweeps,"Quiet");
    printfln("\nGround State Energy = %.10f",energy);
    //printfln("\nGround State Energy pbc = %.10f",energypbc);



    // for(int j = 1; j <= N; ++j)
    //     {
    //     //re-gauge psi to get ready to measure at position j
    //     auto ampo2 = AutoMPO(sites) ;
    //     ampo2 += 1, "Cdag", j, "C", 1;
    //     n_j = toMPO(ampo2);
    //     auto wave1 = inner(psi , n_j , psi) ;
    //     auto wave2 = inner(psi1 , n_j , psi) ;
    //     auto wave3 = inner(psi , n_j , psi1) ;
    //     auto wave4 = inner(psi1 , n_j , psi1) ;
    //     auto b     = inner(psi1,psi) ;
    //     auto majorana = wave1 - b*wave2 - b*wave3 + b*b*wave4 ;
    //     printfln("Site %d occupation: %.5f", j, majorana);
    // }

    //std::vector<double> majorana_values;
    for(int j = 1; j <= N; ++j)
        {
        //re-gauge psi to get ready to measure at position j
        auto ampo3 = AutoMPO(sites) ;
        ampo3 += 1, "C", j, "C", 1;
        G1 = toMPO(ampo3);
        auto ampo4 = AutoMPO(sites) ;
        ampo4 += 1, "C", j, "Cdag", 1;
        G2 = toMPO(ampo4);
        auto ampo5 = AutoMPO(sites) ;
        ampo5 += -1, "C", j, "C", N;
        G3 = toMPO(ampo5);
        auto ampo6 = AutoMPO(sites) ;
        ampo6 += 1, "C", j, "Cdag", N;
        G4 = toMPO(ampo6);
        auto wave1 = inner(psi , G1 , psi) ;
        auto wave2 = inner(psi , G2 , psi) ;
        auto wave3 = inner(psi , G3 , psi) ;
        auto wave4 = inner(psi , G4 , psi) ;
        auto majorana= wave1+ wave2 + wave3 + wave4 ;
        printfln("Site %d occupation: %.5f", j, majorana);
        // printfln("wave1: %.5f",  wave1);
        // printfln("wave2: %.5f",  wave2);
        // printfln("wave3: %.5f",  wave3);
        // printfln("wave4: %.5f",  wave4);
       // majorana_values.push_back(majorana); 
    }
    
    // std::vector<double> result = QML(majorana_values) ;
    // printfln("QML",result);
    return 0;
}


