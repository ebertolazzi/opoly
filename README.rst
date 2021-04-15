OPoly
=====

|File Exchange|

A small header only library for the computation
of orthogonal polynomial.

The supported polynomials are:

- Jacobi
- Legendre
- Chebyshev
- Laguerre
- Hermite

Moreover the class computes the
zeros of the polynomilas and weight and nodes of gaussian quadrature.

Online doc at: http://ebertolazzi.github.io/opoly/

Matlab Interface
----------------

Install the toolbox (download at the
`link <https://github.com/ebertolazzi/opoly/releases>`__).
After the installation compile the Toolbox
in the command windows of Matlab:

.. code-block:: matlab

  > CompileOPolyLib

Usage
~~~~~

Instantiate a polynomial

.. code-block:: matlab

  % kind = 'jacobi', 'legendre', 'chebyshev', 'laguerre', 'hermite'
  kind = 'chebyshev'
  p = OPoly(kind);

special case for Jacobi and Laguerre

.. code-block:: matlab

  p = OPoly('jacobi',alpha,beta);
  p = OPoly('laguerre',alpha);

After instantiation, the polynomial can be evaluated:

.. code-block:: matlab

  n = 10; % polynomial degree
  x = 0:0.1:1;
  y = p.eval( n, x );

evaluate the polynomial with first derivative
and sign variations

.. code-block:: matlab

  [y,yp,s] = p.eval2( n, x );
  % y  = p(x)
  % yp = p'(x)
  % s  = sign variaton of the corresponding Sturm sequence

Its possible to get the complete
Sturm sequence

.. code-block:: matlab

  [p,s] = p.sequence( n, x );
  % p  = p_0(x), p_2(x), ..., p_n(x)
  % s  = sign variaton of the corresponding Sturm sequence

Its possible to compute the zeros:

.. code-block:: matlab

  % epsi = tolerance in the zeros computation
  x = p.zeros( n, epsi );
  % x = vector with the coordinates of the zeros

Finally, compute nodes and weigth of Gauss quadrature
for the interval [-1,1]

.. code-block:: matlab

  % epsi = tolerance in the zeros computation
  [x,w] = p.gauss( n, epsi );
  % x = nodes of the quadrature
  % w = weight of the quadrature

Developer
~~~~~~~~~

| Enrico Bertolazzi
| Dipartimento di Ingegneria Industriale
| Università degli Studi di Trento
| mailto:enrico.bertolazzi@unitn.it (`homepage <www.ing.unitn.it/~bertolaz>`__)


Reference
---------

**Walter Gautschi**,
*Orthogonal Polynomials Computation and Approximation*
Oxford University Press, 2004

.. |File Exchange| image:: https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg
   :target: https://www.mathworks.com/matlabcentral/fileexchange/54481-opoly
