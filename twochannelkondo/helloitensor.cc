#include "itensor/all.h"
#include "itensor/util/print_macro.h"
using namespace itensor ;

int main()
{
//  auto i=Index(3," index i");
//  auto j=Index(4," index j");
//  auto k=Index(5," index k");
//  auto T=randomITensor(i,j,k);
//   T.set(k=3,i=2,j=1, 4.56);
//   T.set(k=1,i=1,j=1, 1);
//   T.set(k=2,i=2,j=1, 2);
// PrintData(T) ;
// PrintData(tags(i)) ;
// auto rho =T * prime(T,i) ;
// Print(rho) ;


// Print(elt(T,k=3,i=2,j=2)) ;
// println("The norm of T is ",norm(T));
//  auto [U,S,V]=svd(T,{i,k});
// //  //Print(norm(U*S*V-T));
// PrintData(U);
// PrintData(S);
// PrintData(V);
// auto m =id("index i") ;
// print(m) ;

int N = 3;
auto sites = SpinHalf(N);
print(sites) ;
// auto state = InitState(sites,"Up");
// auto psi = randomMPS(state);
// auto b=(N+1)/2 ;
// psi.position(b); 

// auto psidag = dag(psi);
// psidag.prime("Link");

// //index linking i to i-1:
// auto li_1 = leftLinkIndex(psi,6);
// // auto rho = prime(psi(6),li_1)*prime(psidag(6),"Site"); 

// auto rho = psi(b) ;
// for(int k=b+1 ; k<=N ;++k)
// {
//  rho *= psi(k) ;
// }
// rho =  prime(rho, "Site")*dag(rho) ;
// // for(int k = 7; k <= 10; ++k)



// //   {
// //   auto lk = rightLinkIndex(psi,k);
// //   rho *= prime(psi(k),lk);
// //   rho *= prime(psidag(k),"Site");
// //   }



// // 移除Link索引


// // Print(rho) ;
// // PrintData(rho(1)) ;

// // auto i = 4;
// // auto j = 10;

// // //'gauge' the MPS to site i
// // //any 'position' between i and j, inclusive, would work here
// // psi.position(i); 

// // auto psidag = dag(psi);
// // psidag.prime("Link");

// // //index linking i to i-1:
// // auto li_1 = leftLinkIndex(psi,i);

// // auto rho = prime(psi(i),li_1)*prime(psidag(i),"Site");
// // for(int k = i+1; k < j; ++k)
// //     {
// //     rho *= psi(k);
// //     rho *= psidag(k);
// //     }
// // //index linking j to j+1:
// // auto lj = rightLinkIndex(psi,j);

// // rho *= prime(psi(j),lj);
// // rho *= prime(psidag(j),"Site");

// //Print(rho) ;

// auto T= inds(rho) ;
// // Print(T) ;




// auto R=T[0] ;
// auto M=T[b] ;
// // Print(R) ;
// // Print(M) ;

//  auto W = M ;
//  M=R ;
//  R=W ;
// // Print(R) ;
// // Print(M) ;



//  rho =replaceInds(rho,{R,M},{M,R}) ;
//  Print(rho) ;

// // auto [U,V] =diagHermitian(rho) ;
// // Print(V) ;
// // Print(U) ;


// // auto u = commonIndex(V,U);
// // Print(u) ;

// // Real Svv = 0.;
// // for(auto n : range1(dim(u)))
// //     {
// //     auto Sn = elt(V,n,n);
// //     Svv += abs(Sn) ;
// //     }
// // Svv= Svv-1 ;
// // Print(Svv) ;

// // auto M= T(1) ;
// // auto N= T(2) ;
// // rho = swapInds(rho,M,N) ;




// // PrintData(T) ;
// // PrintData(M) ;

// // PrintData(psi) ;
// // PrintData(siteIndex(psi,3)) ;
// // PrintData(siteIndex(dag(psi),3)) ;
// // PrintData(tags(siteIndex(psi,3))) ;
// // swapInds(psi(3),)



// // psi.position(2) ;
// // PrintData(psi) ;
// // PrintData(siteIndex(psi,3)) ;
// // PrintData(linkInds(psi)) ;
// // PrintData(leftLinkIndex(psi,3)) ;


// // psi.position(2); 

// // auto psidag = dag(psi);
// // psidag.prime("Link");

// //index linking i to i-1:
// // auto li_1 = leftLinkIndex(psi,2);
// // auto li_2 = siteIndex(psi,2);
// // PrintData(li_1) ;
// // PrintData(li_2) ;
// // auto rho = psi(1)* prime(psidag(1),"Sites") ;
// // for(int k = 2; k <= N; ++k)
// //     {
// //     rho *= psi(k);
// //     rho *= psidag(k);
// //     }
// // Print(rho) ;

  return 0;
}






