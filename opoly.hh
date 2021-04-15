/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  OPOLY                                                                   |
 |                                                                          |
 |  date         : 2000, Aug 27                                             |
 |  release      : 2.0                                                      |
 |  release_date : 1998, Jan 4                                              |
 |  file         : opoly.hh                                                 |
 |  author       : Enrico Bertolazzi                                        |
 |  email        : enrico.bertolazzi@ing.unitn.it                           |
 |                                                                          |
 |  This program is free software; you can redistribute it and/or modify    |
 |  it under the terms of the GNU General Public License as published by    |
 |  the Free Software Foundation; either version 2, or (at your option)     |
 |  any later version.                                                      |
 |                                                                          |
 |  This program is distributed in the hope that it will be useful,         |
 |  but WITHOUT ANY WARRANTY; without even the implied warranty of          |
 |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
 |  GNU General Public License for more details.                            |
 |                                                                          |
 |  You should have received a copy of the GNU General Public License       |
 |  along with this program; if not, write to the Free Software             |
 |  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               |
 |                                                                          |
 |  Copyright (C) 1998                                                      |
 |            ___    ____  ___   _   _        ___    ____  ___   _   _      |
 |           /   \  /     /   \  \  /        /   \  /     /   \  \  /       |
 |          /____/ /__   /____/   \/        /____/ /__   /____/   \/        |
 |         /   \  /     /  \      /        /   \  /     /  \      /         |
 |        /____/ /____ /    \    /        /____/ /____ /    \    /          |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria                                          |
 |      Meccanica e Strutturale               tel: +39-461-882590           |
 |      Universita` degli Studi di Trento     fax: +39-461-882599           |
 |      Via Mesiano 77                                                      |
 |      I-38050 Trento, Italy        email: enrico.bertolazzi@ing.unitn.it  |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# ifndef OPOLY_HH
# define OPOLY_HH

# include <vector>

namespace opoly {

  using std::vector;
  
  # define OPOLY_TYPES_FROM_TYPENAME(T)                            \
  typedef typename T::int_type              int_type;              \
  typedef typename T::index_pointer         index_pointer;         \
  typedef typename T::index_const_pointer   index_const_pointer;   \
  typedef typename T::index_reference       index_reference;       \
  typedef typename T::index_const_reference index_const_reference; \
                                                                   \
  typedef typename T::real_type             real_type;             \
  typedef typename T::pointer               pointer;               \
  typedef typename T::const_pointer         const_pointer;         \
  typedef typename T::reference             reference;             \
  typedef typename T::const_reference       const_reference

  # define OPOLY_TYPES(INT,REAL)                                  \
    typedef INT                int_type;                          \
    typedef int_type *         index_pointer;                     \
    typedef int_type const *   index_const_pointer;               \
    typedef int_type           index_reference;                   \
    typedef int_type const &   index_const_reference;             \
                                                                  \
    typedef REAL               real_type;                         \
    typedef real_type *        pointer;                           \
    typedef real_type const *  const_pointer;                     \
    typedef real_type &        reference;                         \
    typedef real_type const &  const_reference

  template <typename INT, typename REAL>
  class poly {
  public:

    OPOLY_TYPES(INT,REAL);
    poly<INT,REAL> const & operator = (poly<INT,REAL> const &);

  public:
    poly() {};
    virtual ~poly() {};

    real_type
    operator () (int_type n, const_reference x) const {
      real_type p, dp;
      int_type sign_variation = (*this)(n, x, p, dp);
      return p;
    }

    virtual
    int_type operator () (
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const = 0;
  };

  // * * * * * * * * * * * * * * JACOBI * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class Jacobi : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);
    
  private:
    real_type alpha, beta, ab, a2b2;
  
  public:
    Jacobi( const_reference _alpha, const_reference _beta )
    : alpha(_alpha), beta(_beta),
      ab(alpha+beta), a2b2(alpha*alpha-beta*beta)
    { }
    
    int_type
    operator () (
      int_type        n,
      const_reference x,
      reference       p1,
      reference       dp
    ) const override;
    
    Jacobi<INT,REAL> const &
    operator = (Jacobi<INT,REAL> const & P) {
      alpha = P.alpha;
      beta  = P.beta;
      ab    = P.ab;
      a2b2  = P.a2b2;
    }

  };

  template <typename INT, typename REAL>
  inline
  typename Jacobi<INT,REAL>::int_type
  Jacobi<INT,REAL>::operator () (
    int_type        n,
    const_reference x,
    reference       p,
    reference       dp
  ) const {

    real_type p0  = 1, p1  = (alpha-beta)/2+(1+ab/2)*x;
    real_type dp0 = 0, dp1 = 1+ab/2;
        
    real_type oldsign = p1 == 0 ? p0 : p1;
    int_type s       = p1 < 0  ? 1  : 0;

    for ( int_type i=1; i < n; ++i ) {
      real_type a1 = 2*(i+1)*(i+ab+1)*(2*i+ab);
      real_type a2 = (2*i+ab+1)*a2b2;
      real_type a3 = (2*i+ab)*(2*i+ab+1)*(2*i+ab+2);
      real_type a4 = 2*(i+alpha)*(i+beta)*(2*i+ab+2);
      real_type a5 = (a2+x*a3);

      real_type p2  = (          a5 * p1  - a4 * p0 )/a1;
      real_type dp2 = (a3 * p1 + a5 * dp1 - a4 * dp0)/a1;

      p0  = p1; p1  = p2;
      dp0 = dp1; dp1 = dp2;

      if ( oldsign * p1 < 0 ) ++s;
      if ( p1 != 0 ) oldsign = p1;
    }
    p  = p1;
    dp = dp1;
    return s;
  }

  // * * * * * * * * * * * * * * LEGENDRE * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class Legendre : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Legendre() {}
    int_type operator () (
      int_type        n,
      const_reference x,
      reference       p1,
      reference       dp
    ) const override;
  };
  
  template <typename INT, typename REAL>
  inline
  typename Legendre<INT,REAL>::int_type
  Legendre<INT,REAL>::operator () (
    int_type        n,
    const_reference x,
    reference       p,
    reference       dp
  ) const {
    real_type p0  = 1, p1  = x;
    real_type dp0 = 0, dp1 = 1;

    real_type oldsign = p1 == 0 ? p0 : p1;
    int_type s       = p1 < 0  ? 1  : 0;

    for ( int_type i=1; i < n; ++i ) {
      real_type A   = (2*i+1.0)/(i+1.0);
      real_type C   = 1-A;
      real_type p2  =        A*x*p1  + C*p0;
      real_type dp2 = A*p1 + A*x*dp1 + C*dp0;

      p0  = p1; p1  = p2;
      dp0 = dp1; dp1 = dp2;
      
      if ( oldsign * p1 < 0 ) ++s;
      if ( p1 != 0 ) oldsign = p1;
    }
    p  = p1;
    dp = dp1;
    return s;
  }

  // * * * * * * * * * * * * * * CHEBYSHEV * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class Chebyshev : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Chebyshev() {}
    int_type
    operator () (
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const override;
  };
  
  template <typename INT, typename REAL>
  inline
  typename Chebyshev<INT,REAL>::int_type
  Chebyshev<INT,REAL>::operator () (
    int_type        n,
    const_reference x,
    reference       p,
    reference       dp
  ) const {

    real_type p0  = 1, p1  = x;
    real_type dp0 = 0, dp1 = 1;

    real_type oldsign = p1 == 0 ? p0 : p1;
    int_type s       = p1  < 0 ? 1  : 0;

    for ( int_type i=1; i < n; ++i ) {
      real_type p2  =        2*x*p1  - p0;
      real_type dp2 = 2*p1 + 2*x*dp1 - dp0;
      p0  = p1; p1  = p2;
      dp0 = dp1; dp1 = dp2;
      if ( oldsign * p1 < 0 ) ++s;
      if ( p1 != 0 ) oldsign = p1;
    }
    p  = p1;
    dp = dp1;
    return s;
  }

  // * * * * * * * * * * * * * * LAGUERRE * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class Laguerre : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Laguerre() {}
    int_type
    operator () (
      int_type        n,
      const_reference x,
      reference       p1,
      reference       dp
    ) const override;
  };
  
  template <typename INT, typename REAL>
  inline
  typename Laguerre<INT,REAL>::int_type
  Laguerre<INT,REAL>::operator () (
    int_type        n,
    const_reference x,
    reference       p,
    reference       dp
  ) const {

    real_type p0  = 1, p1  = 1-x;
    real_type dp0 = 0, dp1 = -1;
    
    real_type oldsign = p1 == 0 ? p0 : p1;
    int_type  s       = p1  < 0 ? 1  : 0;

    real_type ri = 1;
    for ( int_type i=1; i < n; ++i ) {
      real_type ri1 = ri+1;
      real_type p2  = (      (ri+ri1-x)*p1  - ri*p0 )/ri1;
      real_type dp2 = (-p1 + (ri+ri1-x)*dp1 - ri*dp0)/ri1;
      p0  = p1; p1  = p2;
      dp0 = dp1; dp1 = dp2;
      ri = ri1;
      if ( oldsign * p1 < 0 ) ++s;
      if ( p1 != 0 ) oldsign = p1;
    }
    p  = p1;
    dp = dp1;
    return s;
  }

  // * * * * * * * * * * * * * * HERMITE * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class Hermite : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Hermite() {}
    int_type operator () (
      int_type        n,
      const_reference x,
      reference       p1,
      reference       dp
    ) const override;
  };
  
  template <typename INT, typename REAL>
  inline
  typename Hermite<INT,REAL>::int_type
  Hermite<INT,REAL>::operator () (
    int_type        n,
    const_reference x,
    reference       p,
    reference       dp
  ) const {

    real_type p0  = 1, p1  = 2*x;
    real_type dp0 = 0, dp1 = 2;
        
    real_type oldsign = p1 == 0 ? p0 : p1;
    int_type s       = p1  < 0 ? 1  : 0;

    for ( int_type i=1; i < n; ++i ) {
      real_type p2  = 2*(     x*p1  - i*p0);
      real_type dp2 = 2*(p1 + x*dp1 - i*dp0);
      p0  = p1; p1  = p2;
      dp0 = dp1; dp1 = dp2;
      if ( oldsign * p1 < 0 ) ++s;
      if ( p1 != 0 ) oldsign = p1;
    }
    p  = p1;
    dp = dp1;
    return s;
  }

  // * * * * * * * * * * * * * * ZEROS * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class zeros {
  public:

    OPOLY_TYPES(INT,REAL);

    class interval_type {
    public:      
      real_type a,   // estremo sinistro dell'intervallo
                b;  // estremo desctro dell'intervallo
      int_type  sa,  // variazione di segno in a
                sb; // variazione di segno in b

      interval_type() : a(0), b(0), sa(0), sb(0) {};
   
      interval_type(
        const_reference _a,
        const_reference _b,
        int_type      _sa,
        int_type      _sb)
      : a(_a), b(_b), sa(_sa), sb(_sb) {};

      void
      set(
        const_reference _a,
        const_reference _b,
        int_type        _sa,
        int_type        _sb
      ) {
        a  = _a;
        b  = _b;
        sa = _sa;
        sb = _sb;
      }

      interval_type const &
      operator = (interval_type const & I) {
        a  = I . a;
        b  = I . b;
        sa = I . sa;
        sb = I . sb;
        return *this;
      }

    };
  
    typedef interval_type *       interval_pointer;
    typedef interval_type const * interval_const_pointer;
    typedef interval_type &       interval_reference;
    typedef interval_type const & interval_const_reference;

  private:

    int_type abs_diff(int_type a, int_type b) const
    { return a > b ? a - b : b - a; }

    int_type             N;
    vector<interval_type>  Intervals;
    poly<INT,REAL> const & Poly;
  
    void      sturm_separate( interval_const_reference I );
    void      sturm_interval( interval_reference I ) const;
    real_type sturm_refine( const_reference epsi, interval_const_reference I ) const;
  public:
  
    zeros(poly<INT,REAL> const & pref) : Poly(pref) { }
    ~zeros() { }
    
    // set(POLY const * pP) { ppoly = pP; }
    void eval( int_type n, const_reference epsi, pointer x );
  };

  // * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * *
  
  template <typename INT, typename REAL>
  inline
  void
  zeros<INT,REAL>::sturm_interval(interval_reference I) const {
    // costruisce gli intervalli
    real_type p, pp;
    I.b = 1;
    do {
      I.b  *= 2;
      I.a  = - I.b;
      I.sa = Poly(N, I.a, p, pp);
      I.sb = Poly(N, I.b, p, pp);
    } while ( abs_diff(I.sb, I.sa) != N );
  }

  template <typename INT, typename REAL>
  inline
  typename zeros<INT,REAL>::real_type
  zeros<INT,REAL>::sturm_refine( const_reference epsi, interval_const_reference I ) const {

    real_type p, pp, x;
    int_type count;

    real_type a  = I.a;
    real_type b  = I.b;
    int_type sa = I.sa;

    for ( count = 0; count < 7; ++count ) {
      x = (a + b) / 2;
      int_type sl = Poly(N, x, p, pp);
      if ( sl == sa ) a = x;
      else            b = x;
    } 

    // punto iniziale
    x = (a+b) / 2;
    for ( count = 0; count < 10; ++count ) {
      Poly(N, x, p, pp);
      real_type dx = p/pp;
      x -= dx;
      if ( dx < epsi && dx > -epsi ) break;
    };
    return x;
  }

  // * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * *
  template <typename INT, typename REAL>
  inline
  void
  zeros<INT,REAL>::sturm_separate(interval_const_reference I) {

    switch ( abs_diff(I.sb,I.sa) ) {
    case 0: // se I.sb == I.sa non ci sono radici nell'intervallo [I.a, I.b]
    break;

    case 1: // se I.sb+1 == I.sa c'e` una sola radice nell'intervallo [I.a, I.b]
      Intervals . push_back(I);
    break;
    
    default:
      {
        // suddivido l'intervallo [IN.a,IN.b] in due sottointervalli
        // [IN.a,c] e [c,IN.b] applico la procedura ricorsivamente
        real_type    p, pp, c = (I.a+I.b) / 2;
        int_type    sc = Poly(N, c, p, pp);
        interval_type IA(I.a, c, I.sa, sc);
        interval_type IB(c, I.b, sc, I.sb);
    
        // nell' intervalli [IA.a,IA.b] ci sono na = | IA.sb-IA.sa | radici.
        // verranno isolate negli intervalli I[0], I[1],..., I[na-1]
        sturm_separate(IA);
        // nell' intervalli [IB.a,IB.b] ci sono nb = | IB.sb-IB.sa | radici.
        // verranno isolate negli intervalli I[na+0], I[na+1],..., I[na+nb-1]
        sturm_separate(IB);
      }
    break;
    }
  }

  template <typename INT, typename REAL>
  inline
  void
  zeros<INT,REAL>::eval(int_type n, const_reference epsi, pointer x) {
    N = n;
    Intervals . clear();
    interval_type I;
    sturm_interval(I);
    sturm_separate(I);
    for ( int_type i = 0; i < N; ++i)
      x[i] = sturm_refine(epsi, Intervals[i]);
  }
  
  // * * * * * * * * * * * * * * GaussQ * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL, typename POLY>
  class Gauss_quadrature {
  public:

    OPOLY_TYPES(INT,REAL);
  
  private:
    POLY            lpoly;
    zeros<INT,REAL> zz;
    
  public:
    Gauss_quadrature() : zz(lpoly) {};
    ~Gauss_quadrature() {};

    void
    eval(
      int_type        n,
      const_reference epsi,
      pointer         x,
      pointer         w
    ) {
      zz . eval(n, epsi, x);
      for ( int_type i=0; i < n; ++i) {
        real_type p, pp;
        lpoly(n,x[i],p,pp);
        w[i] = 2/(1-x[i]*x[i])/(pp*pp);
      }
    }
  };
}

# endif
