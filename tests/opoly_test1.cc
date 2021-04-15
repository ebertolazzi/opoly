# include "opoly.hh"

using opoly::Gauss_quadrature;
using opoly::Laguerre;
using opoly::Hermite;
using opoly::Legendre;

typedef unsigned Integer;
typedef double   Real;
typedef opoly::Gauss_quadrature<Integer,Real> GaussQ;

# include <iostream>
# include <iomanip>
# include <string>
# include <cmath>

using namespace std;

Integer const NC   = 15;
Real    const epsi = 0.5/pow(Real(10),Real(NC));

Integer const N = 3;
Real          x[N];
Real          w[N];
GaussQ        gq;

int main() {
  Integer i;
  
  cout.precision(NC);
  gq.eval(N,epsi,x,w);
  
  string banner;
  
  banner = "";
  for ( i = 0; i < NC+5; ++i ) banner += '-';
    banner = "+-----+-" + banner + "-+-" +  banner + "-+";
  
  cout << banner
       << endl
       << "|     | "
       << setw(NC+5) << " coordinates " << " | "
       << setw(NC+5) << " weigths " << " |" 
       << endl
       << banner
       << endl;
  for ( i = 0; i < N; ++i )
    cout << "| "
         << setw(3) << i << " | "
         << setw(NC+5) << x[i] << " | "
         << setw(NC+5) << w[i] << " |" 
         << endl;
  cout << banner << endl;
}

