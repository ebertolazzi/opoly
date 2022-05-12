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

#ifndef OPOLY_HH
#define OPOLY_HH

#include <vector>
#include <cmath>
#include <string>

namespace opoly {

  using std::vector;
  using std::pow;
  using std::abs;
  using std::string;

  # define OPOLY_TYPES_FROM_TYPENAME(T)                        \
  typedef typename T::int_type            int_type;            \
  typedef typename T::int_pointer         int_pointer;         \
  typedef typename T::int_const_pointer   int_const_pointer;   \
  typedef typename T::int_reference       int_reference;       \
  typedef typename T::int_const_reference int_const_reference; \
                                                               \
  typedef typename T::real_type           real_type;           \
  typedef typename T::pointer             pointer;             \
  typedef typename T::const_pointer       const_pointer;       \
  typedef typename T::reference           reference;           \
  typedef typename T::const_reference     const_reference

  # define OPOLY_TYPES(INT,REAL)                              \
    typedef INT                int_type;                      \
    typedef int_type *         int_pointer;                   \
    typedef int_type const *   int_const_pointer;             \
    typedef int_type           int_reference;                 \
    typedef int_type const &   int_const_reference;           \
                                                              \
    typedef REAL               real_type;                     \
    typedef real_type *        pointer;                       \
    typedef real_type const *  const_pointer;                 \
    typedef real_type &        reference;                     \
    typedef real_type const &  const_reference

  //!
  //! Base class defining an orthogonal polynomial
  //!
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //!
  template <typename INT, typename REAL>
  class poly {

    poly<INT,REAL> const & operator = (poly<INT,REAL> const &) = delete;

  public:

    OPOLY_TYPES(INT,REAL);

    poly() {};

    virtual ~poly() {};

    virtual string type() const = 0;

    //!
    //! Evaluate the polynomial
    //!
    //! \param[in] n the degree op the polynomial
    //! \param[in] x the value at which the polynomial is evaluated
    //! \return the value \f$ p(x) \f$
    //!
    virtual
    real_type
    operator () ( int_type n, const_reference x ) const = 0;

    //!
    //! Evaluate the sign variation of the corresponding sturm sequence.
    //!
    //! \param[in] n   the degree op the polynomial
    //! \param[in] x   the value at which the polynomial is evaluated
    //! \return    the number of sign variations
    //!
    virtual
    int_type
    svar( int_type n, const_reference x ) const = 0;

    //!
    //! Evaluate the polynomial and the sign variation of the
    //! corresponding sturm sequence.
    //!
    //! \param[in]  n  the degree op the polynomial
    //! \param[in]  x  the value at which the polynomial is evaluated
    //! \param[out] p  the value \f$ p(x) \f$
    //! \return the number of sign variations
    //!
    virtual
    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p
    ) const = 0;

    //!
    //! Evaluate the polynomial and its derivative.
    //! Moreover return the sign variation of the
    //! corresponding sturm sequence.
    //!
    //! \param[in]  n  the degree op the polynomial
    //! \param[in]  x  the value at which the polynomial is evaluated
    //! \param[out] p  the value \f$ p(x) \f$
    //! \param[out] dp the value \f$ p'(x) \f$
    //! \return the number of sign variations
    //!
    virtual
    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const = 0;

    //!
    //! Evaluate the Sturm sequence
    //!
    //! \param[in]  n the degree of the polynomial
    //! \param[in]  x the value at which the polynomial is evaluated
    //! \param[out] pvec the sequence \f$ p_0(x) \f$,  \f$ p_1(x) \f$, ...,  \f$ p_n(x) \f$ of the 3 term recurrence
    //! \return the number of sign variations
    //!
    virtual
    int_type
    sequence( int_type n, const_reference x, pointer pvec ) const = 0;

    //!
    //! Weight function of the orthogonal polynomial
    //! \param[in]  x  the value at which the weight is evaluated
    //!
    virtual
    real_type
    weight( const_reference x ) const = 0;

  };

  // * * * * * * * * * * * * * * JACOBI * * * * * * * * * * * * * * * * * * *

  //!
  //! Class defining the Jacobi orthogonal polynomial
  //!
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //!
  template <typename INT, typename REAL>
  class Jacobi : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);
    
  private:
    real_type alpha, beta, ab, a2b2;
    Jacobi() = delete;

  public:

    Jacobi( const_reference _alpha, const_reference _beta )
    : alpha(_alpha)
    , beta(_beta)
    , ab(_alpha+_beta)
    , a2b2(_alpha*_alpha-_beta*_beta)
    { }

    Jacobi( Jacobi<INT,REAL> const * ptr )
    : alpha(ptr->alpha)
    , beta(ptr->beta)
    , ab(ptr->ab)
    , a2b2(ptr->a2b2)
    { }

    string type() const override { return "jacobi"; }

    //!
    //! Setup the orthogonal polynomial of Jacobi
    //!
    void
    setup( const_reference _alpha, const_reference _beta ) {
      alpha = _alpha;
      beta  = _beta;
      ab    = _alpha+_beta;
      a2b2  = _alpha*_alpha -_beta*_beta;
    }

    real_type
    operator () ( int_type n, const_reference x ) const override {

      real_type p0 = 1;
      real_type p1 = (alpha-beta)/2+(1+ab/2)*x;

      for ( int_type i=1; i < n; ++i ) {
        real_type t1 = 2*i+ab;
        real_type a1 = 2*(i+1)*(i+ab+1)*t1;
        real_type a2 = (t1+1)*a2b2;
        real_type a3 = t1*(t1+1)*(t1+2);
        real_type a4 = 2*(i+alpha)*(i+beta)*(t1+2);
        real_type a5 = (a2+x*a3);
        real_type p2 = ( a5 * p1 - a4 * p0 )/a1;
        p0 = p1; p1 = p2;
      }
      return p1;
    }

    int_type
    svar( int_type n, const_reference x ) const override {

      real_type p0 = 1;
      real_type p1 = (alpha-beta)/2+(1+ab/2)*x;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1 < 0  ? 1  : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type t1 = 2*i+ab;
        real_type a1 = 2*(i+1)*(i+ab+1)*t1;
        real_type a2 = (t1+1)*a2b2;
        real_type a3 = t1*(t1+1)*(t1+2);
        real_type a4 = 2*(i+alpha)*(i+beta)*(t1+2);
        real_type a5 = (a2+x*a3);

        real_type p2  = ( a5 * p1 - a4 * p0 )/a1;
        p0 = p1; p1 = p2;

        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      return s;
    }

    int_type
    eval( int_type n, const_reference x, reference p ) const override {
      real_type p0 = 1;
      real_type p1 = (alpha-beta)/2+(1+ab/2)*x;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1 < 0  ? 1  : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type t1 = 2*i+ab;
        real_type a1 = 2*(i+1)*(i+ab+1)*t1;
        real_type a2 = (t1+1)*a2b2;
        real_type a3 = t1*(t1+1)*(t1+2);
        real_type a4 = 2*(i+alpha)*(i+beta)*(t1+2);
        real_type a5 = (a2+x*a3);

        real_type p2 = ( a5 * p1 - a4 * p0 )/a1;
        p0 = p1; p1 = p2;

        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p  = p1;
      return s;
    }

    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const override {

      real_type p0  = 1, p1  = (alpha-beta)/2+(1+ab/2)*x;
      real_type dp0 = 0, dp1 = 1+ab/2;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1 < 0  ? 1  : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type t1 = 2*i+ab;
        real_type a1 = 2*(i+1)*(i+ab+1)*t1;
        real_type a2 = (t1+1)*a2b2;
        real_type a3 = t1*(t1+1)*(t1+2);
        real_type a4 = 2*(i+alpha)*(i+beta)*(t1+2);
        real_type a5 = (a2+x*a3);

        real_type p2  = (           a5 * p1  - a4 * p0  )/a1;
        real_type dp2 = ( a3 * p1 + a5 * dp1 - a4 * dp0 )/a1;

        p0  = p1;  p1  = p2;
        dp0 = dp1; dp1 = dp2;

        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p  = p1;
      dp = dp1;
      return s;
    }

    int_type
    sequence( int_type n, const_reference x, pointer p ) const override {
      p[0] = 1;
      p[1] = (alpha-beta)/2+(1+ab/2)*x;

      real_type oldsign = p[1] == 0 ? p[0] : p[1];
      int_type  s       = p[1] < 0  ? 1  : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type t1 = 2*i+ab;
        real_type a1 = 2*(i+1)*(i+ab+1)*t1;
        real_type a2 = (t1+1)*a2b2;
        real_type a3 = t1*(t1+1)*(t1+2);
        real_type a4 = 2*(i+alpha)*(i+beta)*(t1+2);
        real_type a5 = (a2+x*a3);

        real_type new_p = ( a5 * p[i] - a4 * p[i-1] )/a1;

        if ( oldsign * new_p < 0 ) ++s;
        if ( new_p != 0 ) oldsign = new_p;
        p[i+1] = new_p;
      }
      return s;
    }

    Jacobi<INT,REAL> const &
    operator = ( Jacobi<INT,REAL> const & P ) {
      alpha = P.alpha;
      beta  = P.beta;
      ab    = P.ab;
      a2b2  = P.a2b2;
    }

    real_type
    weight( const_reference x ) const override {
      return pow(1-x,alpha)*pow(1+x,beta);
    }
  };

  // * * * * * * * * * * * * * * LEGENDRE * * * * * * * * * * * * * * * * * * *

  //!
  //! Class defining the Legendre orthogonal polynomial
  //!
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //!

  template <typename INT, typename REAL>
  class Legendre : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Legendre() {}

    Legendre( Legendre<INT,REAL> const * ) {}

    string type() const override { return "legendre"; }

    real_type
    operator () ( int_type n, const_reference x ) const override {
      real_type p0 = 1;
      real_type p1 = x;
      for ( int_type i=1; i < n; ++i ) {
        real_type A  = (2*i+1.0)/(i+1.0);
        real_type C  = 1-A;
        real_type p2 = A*x*p1 + C*p0;
        p0 = p1;
        p1 = p2;
      }
      return p1;
    }

    int_type
    svar( int_type n, const_reference x ) const override {
      real_type p0 = 1, p1 = x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1 < 0  ? 1  : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type A  = (2*i+1.0)/(i+1.0);
        real_type C  = 1-A;
        real_type p2 = A*x*p1 + C*p0;
        p0 = p1; p1 = p2;

        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      return s;
    }

    int_type
    eval( int_type n, const_reference x, reference p ) const override {
      real_type p0 = 1, p1 = x;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1 < 0  ? 1  : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type A  = (2*i+1.0)/(i+1.0);
        real_type C  = 1-A;
        real_type p2 = A*x*p1 + C*p0;
        p0 = p1; p1 = p2;

        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p = p1;
      return s;
    }

    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const override {
      real_type p0  = 1, p1  = x;
      real_type dp0 = 0, dp1 = 1;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1 < 0  ? 1  : 0;

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

    int_type
    sequence( int_type n, const_reference x, pointer p ) const override {
      p[0] = 1;
      p[1] = x;

      real_type oldsign = p[1] == 0 ? p[0] : p[1];
      int_type  s       = p[1] < 0  ? 1    : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type A = (2*i+1.0)/(i+1.0);
        real_type C = 1-A;
        real_type new_p = A*x*p[i]  + C*p[i-1];

        if ( oldsign * new_p < 0 ) ++s;
        if ( new_p != 0 ) oldsign = new_p;
        p[i+1] = new_p;
      }
      return s;
    }

    real_type
    weight( const_reference ) const override {
      return 1;
    }

  };

  // * * * * * * * * * * * * * * CHEBYSHEV * * * * * * * * * * * * * * * * * * *

  //!
  //! Class defining the Chebyshev orthogonal polynomial
  //!
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //!
  template <typename INT, typename REAL>
  class Chebyshev : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Chebyshev() {}

    Chebyshev( Chebyshev<INT,REAL> const * ) {}

    string type() const override { return "chebyshev"; }

    real_type
    operator () ( int_type n, const_reference x ) const override {
      real_type p0 = 1, p1 = x;
      for ( int_type i=1; i < n; ++i ) {
        real_type p2 = 2*x*p1 - p0;
        p0 = p1; p1 = p2;
      }
      return p1;
    }

    int_type
    svar( int_type n, const_reference x ) const override {
      real_type p0 = 1, p1 = x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;
      for ( int_type i=1; i < n; ++i ) {
        real_type p2 = 2*x*p1 - p0;
        p0 = p1; p1 = p2;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      return s;
    }

    int_type
    eval( int_type n, const_reference x, reference p ) const override {
      real_type p0 = 1, p1 = x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;
      for ( int_type i=1; i < n; ++i ) {
        real_type p2 = 2*x*p1 - p0;
        p0 = p1; p1 = p2;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p = p1;
      return s;
    }

    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const override {

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

    int_type
    sequence( int_type n, const_reference x, pointer p ) const override {

      p[0] = 1;
      p[1] = x;

      real_type oldsign = p[1] == 0 ? p[0] : p[1];
      int_type  s       = p[1]  < 0 ? 1    : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type new_p = 2*x*p[i] - p[i-1];
        if ( oldsign * new_p < 0 ) ++s;
        if ( new_p != 0 ) oldsign = new_p;
        p[i+1] = new_p;
      }
      return s;
    }

    real_type
    weight( const_reference x ) const override {
      return 1/sqrt(1-x*x);
    }
  };

  // * * * * * * * * * * * * * * LAGUERRE * * * * * * * * * * * * * * * * * * *

  //!
  //! Class defining the Laguerre orthogonal polynomial
  //!
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //!

  template <typename INT, typename REAL>
  class Laguerre : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

  private:
    real_type alpha;

  public:

    Laguerre( real_type _alpha ) : alpha(_alpha) { }

    Laguerre( Laguerre<INT,REAL> const * ptr ) : alpha(ptr->alpha) { }

    string type() const override { return "laguerre"; }

    real_type
    operator () ( int_type n, const_reference x ) const override {
      real_type p0 = 1, p1 = 1-x;
      real_type ri = 1;
      for ( int_type i=1; i < n; ++i ) {
        real_type ri1 = ri+1;
        real_type p2  = ( (ri+ri1+alpha-x)*p1 - (ri+alpha)*p0 )/ri1;
        p0 = p1;
        p1 = p2;
        ri = ri1;
      }
      return p1;
    }

    int_type
    svar( int_type n, const_reference x ) const override {

      real_type p0 = 1, p1 = 1-x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;

      real_type ri = 1;
      for ( int_type i=1; i < n; ++i ) {
        real_type ri1 = ri+1;
        real_type p2  = ( (ri+ri1+alpha-x)*p1 - (ri+alpha)*p0 )/ri1;
        p0 = p1; p1 = p2;
        ri = ri1;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      return s;
    }

    int_type
    eval( int_type n, const_reference x, reference p ) const override {
      real_type p0 = 1, p1 = 1-x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;
      real_type ri = 1;
      for ( int_type i=1; i < n; ++i ) {
        real_type ri1 = ri+1;
        real_type p2 = ( (ri+ri1+alpha-x)*p1 - (ri+alpha)*p0 )/ri1;
        p0 = p1; p1 = p2;
        ri = ri1;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p = p1;
      return s;
    }

    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const override {

      real_type p0  = 1, p1  = 1-x;
      real_type dp0 = 0, dp1 = -1;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;

      real_type ri = 1;
      for ( int_type i=1; i < n; ++i ) {
        real_type ri1 = ri+1;
        real_type p2  = (      (ri+ri1+alpha-x)*p1  - (ri+alpha)*p0 )/ri1;
        real_type dp2 = (-p1 + (ri+ri1+alpha-x)*dp1 - (ri+alpha)*dp0)/ri1;
        p0  = p1;  p1  = p2;
        dp0 = dp1; dp1 = dp2;
        ri = ri1;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p  = p1;
      dp = dp1;
      return s;
    }

    int_type
    sequence( int_type n, const_reference x, pointer p ) const override {

      p[0] = 1;
      p[1] = 1-x;

      real_type oldsign = p[1] == 0 ? p[0] : p[1];
      int_type  s       = p[1]  < 0 ? 1  : 0;

      real_type ri = 1;
      for ( int_type i=1; i < n; ++i ) {
        real_type ri1 = ri+1;
        real_type new_p = ( (ri+ri1+alpha-x)*p[i] - (ri+alpha)*p[i-1] )/ri1;
        ri = ri1;
        if ( oldsign * new_p < 0 ) ++s;
        if ( new_p != 0 ) oldsign = new_p;
        p[i+1] = new_p;
      }
      return s;
    }

    real_type
    weight( const_reference x ) const override {
      return exp(-x)*pow(x,alpha);
    }
  };

  // * * * * * * * * * * * * * * HERMITE * * * * * * * * * * * * * * * * * * *

  //!
  //! Class defining the Hermite orthogonal polynomial
  //!
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //!

  template <typename INT, typename REAL>
  class Hermite : public poly<INT,REAL> {
  public:
    OPOLY_TYPES(INT,REAL);

    Hermite() {}

    Hermite( Hermite<INT,REAL> const *) { }

    string type() const override { return "hermite"; }

    real_type
    operator () ( int_type n, const_reference x ) const override {
      real_type p0 = 1, p1 = 2*x;
      for ( int_type i=1; i < n; ++i ) {
        real_type p2 = 2*( x*p1 - i*p0);
        p0 = p1; p1 = p2;
      }
      return p1;
    }

    int_type
    svar( int_type n, const_reference x ) const override {
      real_type p0 = 1, p1 = 2*x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;
      for ( int_type i=1; i < n; ++i ) {
        real_type p2  = 2*( x*p1 - i*p0 );
        p0  = p1; p1  = p2;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      return s;
    }

    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p
    ) const override {
      real_type p0 = 1, p1 = 2*x;
      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;
      for ( int_type i=1; i < n; ++i ) {
        real_type p2 = 2*( x*p1 - i*p0 );
        p0 = p1; p1 = p2;
        if ( oldsign * p1 < 0 ) ++s;
        if ( p1 != 0 ) oldsign = p1;
      }
      p = p1;
      return s;
    }

    int_type
    eval(
      int_type        n,
      const_reference x,
      reference       p,
      reference       dp
    ) const override {

      real_type p0  = 1, p1  = 2*x;
      real_type dp0 = 0, dp1 = 2;

      real_type oldsign = p1 == 0 ? p0 : p1;
      int_type  s       = p1  < 0 ? 1  : 0;

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

    int_type
    sequence( int_type n, const_reference x, pointer p ) const override {

      p[0] = 1;
      p[1] = 2*x;

      real_type oldsign = p[1] == 0 ? p[0] : p[1];
      int_type  s       = p[1]  < 0 ? 1    : 0;

      for ( int_type i=1; i < n; ++i ) {
        real_type new_p = 2*( x*p[i] - i*p[i-1]);
        if ( oldsign * new_p < 0 ) ++s;
        if ( new_p != 0 ) oldsign = new_p;
      }
      return s;
    }

    real_type
    weight( const_reference x ) const override {
      return exp(-x*x);
    }
  };

  // * * * * * * * * * * * * * * ZEROS * * * * * * * * * * * * * * * * * * *

  //! 
  //! Class computeing the zeros of an orthogonal polynomial
  //! 
  //! - INT  is a integer type, can be also unsigned
  //! - REAL is a floating point type, can be float or a high precision number class
  //! 

  template <typename INT, typename REAL>
  class zeros {
  private:
    zeros() = delete;
    zeros<INT,REAL> const & operator = ( zeros<INT,REAL> const & ) = delete;
  public:

    OPOLY_TYPES(INT,REAL);

    class interval_type {
    public:      
      real_type a;   // estremo sinistro dell'intervallo
      real_type b;   // estremo desctro dell'intervallo
      int_type  sa;  // variazione di segno in a
      int_type  sb;  // variazione di segno in b

      interval_type() : a(0), b(0), sa(0), sb(0) {};
   
      interval_type(
        const_reference _a,
        const_reference _b,
        int_type        _sa,
        int_type        _sb)
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

    int_type               N;
    vector<interval_type>  Intervals;
    poly<INT,REAL> const & Poly;

    void
    sturm_separate( interval_const_reference I ) {

      switch ( abs_diff(I.sb,I.sa) ) {
      case 0: // se I.sb == I.sa non ci sono radici nell'intervallo [I.a, I.b]
      break;

      case 1: // se I.sb+1 == I.sa c'e` una sola radice nell'intervallo [I.a, I.b]
        Intervals.push_back(I);
      break;

      default:
        {
          // suddivido l'intervallo [IN.a,IN.b] in due sottointervalli
          // [IN.a,c] e [c,IN.b] applico la procedura ricorsivamente
          real_type c  = (I.a+I.b) / 2;
          int_type  sc = Poly.svar( N, c );
          interval_type IA( I.a, c, I.sa, sc );
          interval_type IB( c, I.b, sc, I.sb );

          // nell' intervalli [IA.a,IA.b] ci sono na = | IA.sb-IA.sa | radici.
          // verranno isolate negli intervalli I[0], I[1],..., I[na-1]
          sturm_separate( IA );
          // nell' intervalli [IB.a,IB.b] ci sono nb = | IB.sb-IB.sa | radici.
          // verranno isolate negli intervalli I[na+0], I[na+1],..., I[na+nb-1]
          sturm_separate( IB );
        }
      break;
      }
    }

    void
    sturm_interval( interval_reference I ) const {
      // costruisce gli intervalli
      I.b = 1;
      do {
        I.b  *= 2;
        I.a  = - I.b;
        I.sa = Poly.svar(N, I.a);
        I.sb = Poly.svar(N, I.b);
      } while ( abs_diff(I.sb, I.sa) != N );
    }

    real_type
    sturm_refine( const_reference epsi, interval_const_reference I ) const {

      real_type x;
      int_type  count;

      real_type a = I.a;
      real_type b = I.b;
      int_type sa = I.sa;

      for ( count = 0; count < 7; ++count ) {
        x = (a + b) / 2;
        int_type sl = Poly.svar(N, x);
        if ( sl == sa ) a = x;
        else            b = x;
      }

      // punto iniziale
      x = (a + b) / 2;
      for ( count = 0; count < 10; ++count ) {
        real_type p, pp;
        Poly.eval(N, x, p, pp);
        real_type dx = p/pp;
        x -= dx;
        if ( dx < epsi && dx > -epsi ) break;
      };
      return x;
    }

  public:

    zeros( poly<INT,REAL> const & pref ) : Poly(pref) { }
    ~zeros() { }

    void
    eval( int_type n, const_reference epsi, pointer x ) {
      N = n;
      Intervals.clear();
      interval_type I;
      sturm_interval(I);
      sturm_separate(I);
      for ( int_type i = 0; i < N; ++i)
        x[i] = sturm_refine( epsi, Intervals[i]) ;
    }

  };

  // * * * * * * * * * * * * * * GaussQ * * * * * * * * * * * * * * * * * * *

  template <typename INT, typename REAL>
  class Gauss_quadrature {
  public:

    OPOLY_TYPES(INT,REAL);
  
  private:
    Legendre<INT,REAL> lpoly;
    zeros<INT,REAL>    zz;
    
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
      zz.eval( n, epsi, x );
      for ( int_type i=0; i < n; ++i ) {
        real_type p, pp;
        lpoly.eval(n,x[i],p,pp);
        w[i] = 2/(1-x[i]*x[i])/(pp*pp);
      }
    }
  };

}

# endif
