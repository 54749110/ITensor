#include "itensor/all.h"
#include "itensor/util/print_macro.h"
using namespace itensor;
using channelkondo= twochannelkondo<SpinHalfSite,ElectronSite>;

int 
main()
    {
    int N = 100;
    channelkondo sites;
    MPO H;
    MPS psi0;
    //
    // Initialize the site degrees of freedom
    // Setting "ConserveQNs=",true makes the indices
    // carry Sz quantum numbers and will lead to 
    // block-sparse MPO and MPS tensors
    //
    sites = channelkondo(N); //make a chain of N spin 1/2's
    // auto sites = SpinOne(N); //make a chain of N spin 1's
  //  double  k[10]={0,0.1,0.2,0.3,0.5,1,2,3,5,10};
  //  for (auto m=-1;m<9;m+=1)
  //  {
  //  auto V=k[m] ;

   auto  V=1; 
    //
    auto ampo = AutoMPO(sites);
    for(int j=1; j<=N-2 ; ++j)
        {
        ampo += 0.5*V,"S+",1,"Cdagup",j+1,"Cdn",j+1;
        ampo += 0.5*V,"S+",1,"Cdagdn",j+1,"Cup",j+1;
	ampo += 0.5*V,"S+",N,"Cdagup",j+1,"Cdn",j+1;
	ampo += 0.5*V,"S+",N,"Cdagdn",j+1,"Cup",j+1;
	ampo +=-0.5*V,"S+",1,"Cdagup",j+1,"Cdn",j+1;
        ampo +=-0.5*V,"S+",1,"Cdagdn",j+1,"Cup",j+1;
        ampo +=-0.5*V,"S+",N,"Cdagup",j+1,"Cdn",j+1;
        ampo +=-0.5*V,"S+",N,"Cdagdn",j+1,"Cup",j+1;
	ampo += 0.5*V,"S-",1,"Cdagup",j+1,"Cdn",j+1;
        ampo += 0.5*V,"S-",1,"Cdagdn",j+1,"Cup",j+1;
        ampo += 0.5*V,"S-",N,"Cdagup",j+1,"Cdn",j+1;
        ampo += 0.5*V,"S-",N,"Cdagdn",j+1,"Cup",j+1;
        ampo += 0.5*V,"S-",1,"Cdagup",j+1,"Cdn",j+1;
        ampo += 0.5*V,"S-",1,"Cdagdn",j+1,"Cup",j+1;
        ampo += 0.5*V,"S-",N,"Cdagup",j+1,"Cdn",j+1;
        ampo += 0.5*V,"S-",N,"Cdagdn",j+1,"Cup",j+1;
	ampo += V,"Sz",1,"Cdagup",j+1,"Cup",j+1;
        ampo +=-V,"Sz",1,"Cdagdn",j+1,"Cdn",j+1;
        ampo += V,"Sz",N,"Cdagup",j+1,"Cup",j+1;
        ampo +=-V,"Sz",N,"Cdagdn",j+1,"Cdn",j+1;
	ampo +=-0.01,"Cdagdn",j+1,"Cdn",j+1;
        ampo +=-0.01,"Cdagup",j+1,"Cup",j+1;
       	}
    for(int k=1; k<=N-3; ++k)
      {
        ampo += 1,"Cdagup",k+1,"Cdn",k+2;
	ampo += 1,"Cdagdn",k+1,"Cup",k+2;
      }
    // ampo += -0.01,"Cdagup",N-1,"Cup",N-1;
    // ampo += -0.01,"Cdagdn",N-1,"Cdn",N-1;
    H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(auto i : range1(N))
        {
        if(i%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    psi0 = MPS(state);

    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    println("Initial energy = %.5f", inner(psi0,H,psi0) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    println("\nGround State Energy = %.10f",energy);
    println("\nUsing inner = %.10f", inner(psi,H,psi) );
   //  }
    return 0;
    }
