#include "itensor/all.h"
#include "itensor/util/print_macro.h"
using namespace itensor ;

int main()
{
auto i=Index(3," index i");
auto j=Index(4," index j");
auto k=Index(5," index k");
auto T=randomITensor(i,j,k);
auto [U,S,V]=svd(T,{i,k});
Print(norm(U*S*V-T));

return 0;
}






