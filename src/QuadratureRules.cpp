#include "QuadratureRules.hpp"

//----------------------------------- core implementation 

struct GaussianQuadratureRules
{
  GaussianQuadratureRules (void) {}

  // Everything in this struct is slightly adapted from
  //   https://people.sc.fsu.edu/~jburkardt/cpp_src/gen_laguerre_rule/gen_laguerre_rule.cpp
  // under the terms of the GNU GPL. Thanks, Prof. Burkardt.

  void cdgqf (int nt, int kind, double alpha, double beta, double t[], double wts[])
  {
    //  Purpose:
    //    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
    //
    //  Discussion:
    //    This routine computes all the knots and weights of a Gauss quadrature
    //    formula with a classical weight function with default values for A and B,
    //    and only simple knots.
    //
    //    There are no moments checks and no printing is done.
    //
    //    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    08 January 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //  Parameters:
    //    Input, int NT, the number of knots.
    //    Input, int KIND, the rule.
    //    1, Legendre,             (a,b)       1.0
    //    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
    //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    //    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    //    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    //    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
    //    Input, double ALPHA, the value of Alpha, if needed.
    //    Input, double BETA, the value of Beta, if needed.
    //    Output, double T[NT], the knots.
    //    Output, double WTS[NT], the weights.

    std::cout << "cdgqf\n" << std::endl;//DEBUG

    double *aj;
    double *bj;
    double zemu;

    parchk ( kind, 2 * nt, alpha, beta );

    // Get the Jacobi matrix and zero-th moment.
    aj = new double[nt];
    bj = new double[nt];

    zemu = class_matrix ( kind, nt, alpha, beta, aj, bj );

    // Compute the knots and weights.
    sgqf ( nt, aj, bj, zemu, t, wts );

    delete [] aj;
    delete [] bj;

    std::cout << "leave cdgqf\n" << std::endl;//DEBUG

    return;
  }

  void cgqf (int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[])
  {
    //  Purpose:
    //    CGQF computes knots and weights of a Gauss quadrature formula.
    //
    //  Discussion:
    //    The user may specify the interval (A,B).
    //
    //    Only simple knots are produced.
    //
    //    Use routine EIQFS to evaluate this quadrature formula.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    16 February 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //  Parameters:
    //    Input, int NT, the number of knots.
    //    Input, int KIND, the rule.
    //    1, Legendre,             (a,b)       1.0
    //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
    //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
    //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
    //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
    //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
    //    Input, double ALPHA, the value of Alpha, if needed.
    //    Input, double BETA, the value of Beta, if needed.
    //    Input, double A, B, the interval endpoints, or other parameters.
    //    Output, double T[NT], the knots.
    //    Output, double WTS[NT], the weights.

    std::cout << "cgqf\n" << std::endl;//DEBUG

    int i;
    int *mlt;
    int *ndx;

    // Compute the Gauss quadrature formula for default values of A and B.
    cdgqf ( nt, kind, alpha, beta, t, wts );

    // Prepare to scale the quadrature formula to other weight function with valid A and B.
    mlt = new int[nt];
    for ( i = 0; i < nt; i++ )
    {
      mlt[i] = 1;
    }
    ndx = new int[nt];
    for ( i = 0; i < nt; i++ )
    {
      ndx[i] = i + 1;
    }
    scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b );

    delete [] mlt;
    delete [] ndx;

    std::cout << "leave cgqf\n" << std::endl;//DEBUG

    return;
  }

  double class_matrix (int kind, int m, double alpha, double beta, double aj[], double bj[])
  {
    //  Purpose:
    //    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
    //
    //  Discussion:
    //    This routine computes the diagonal AJ and sub-diagonal BJ
    //    elements of the order M tridiagonal symmetric Jacobi matrix
    //    associated with the polynomials orthogonal with respect to
    //    the weight function specified by KIND.
    //
    //    For weight functions 1-7, M elements are defined in BJ even
    //    though only M-1 are needed.  For weight function 8, BJ(M) is
    //    set to zero.
    //
    //    The zero-th moment of the weight function is returned in ZEMU.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    08 January 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //  Parameters:
    //    Input, int KIND, the rule.
    //    1, Legendre,             (a,b)       1.0
    //    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
    //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    //    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    //    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    //    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
    //    Input, int M, the order of the Jacobi matrix.
    //    Input, double ALPHA, the value of Alpha, if needed.
    //    Input, double BETA, the value of Beta, if needed.
    //    Output, double AJ[M], BJ[M], the diagonal and subdiagonal of the Jacobi matrix.
    //    Output, double CLASS_MATRIX, the zero-th moment.

    std::cout << "class matrix\n" << std::endl;//DEBUG

    double a2b2;
    double ab;
    double aba;
    double abi;
    double abj;
    double abti;
    double apone;
    int i;
    double pi = 3.14159265358979323846264338327950;
    double temp;
    double temp2;
    double zemu;

    temp = r8_epsilon ( );

    parchk ( kind, 2 * m - 1, alpha, beta );

    temp2 = 0.5;

    if ( 500.0 * temp < std::fabs ( std::pow ( std::tgamma ( temp2 ), 2 ) - pi ) )
      Rcpp::stop("CLASS_MATRIX: Fatal error.\n\tGamma function does not match machine parameters.");

    if ( kind == 1 )
    {
      ab = 0.0;

      zemu = 2.0 / ( ab + 1.0 );

      for ( i = 0; i < m; i++ )
      {
        aj[i] = 0.0;
      }

      for ( i = 1; i <= m; i++ )
      {
        abi = i + ab * ( i % 2 );
        abj = 2 * i + ab;
        bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
      }
    }
    else if ( kind == 2 )
    {
      zemu = pi;

      for ( i = 0; i < m; i++ )
      {
        aj[i] = 0.0;
      }

      bj[0] =  sqrt ( 0.5 );
      for ( i = 1; i < m; i++ )
      {
        bj[i] = 0.5;
      }
    }
    else if ( kind == 3 )
    {
      ab = alpha * 2.0;
      zemu = std::pow ( 2.0, ab + 1.0 ) * std::pow ( std::tgamma ( alpha + 1.0 ), 2 )
        / std::tgamma ( ab + 2.0 );

      for ( i = 0; i < m; i++ )
      {
        aj[i] = 0.0;
      }

      bj[0] = sqrt ( 1.0 / ( 2.0 * alpha + 3.0 ) );
      for ( i = 2; i <= m; i++ )
      {
        bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.0 * pow ( i + alpha, 2 ) - 1.0 ) );
      }
    }
    else if ( kind == 4 )
    {
      ab = alpha + beta;
      abi = 2.0 + ab;
      zemu = std::pow ( 2.0, ab + 1.0 ) * std::tgamma ( alpha + 1.0 ) 
        * std::tgamma ( beta + 1.0 ) / std::tgamma ( abi );
      aj[0] = ( beta - alpha ) / abi;
      bj[0] = sqrt ( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
          / ( ( abi + 1.0 ) * abi * abi ) );
      a2b2 = beta * beta - alpha * alpha;

      for ( i = 2; i <= m; i++ )
      {
        abi = 2.0 * i + ab;
        aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
        abi = abi * abi;
        bj[i-1] = sqrt ( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) 
            / ( ( abi - 1.0 ) * abi ) );
      }
    }
    else if ( kind == 5 )
    {
      zemu = std::tgamma ( alpha + 1.0 );

      for ( i = 1; i <= m; i++ )
      {
        aj[i-1] = 2.0 * i - 1.0 + alpha;
        bj[i-1] = sqrt ( i * ( i + alpha ) );
      }
    }
    else if ( kind == 6 )
    {
      zemu = std::tgamma ( ( alpha + 1.0 ) / 2.0 );

      for ( i = 0; i < m; i++ )
      {
        aj[i] = 0.0;
      }

      for ( i = 1; i <= m; i++ )
      {
        bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.0 );
      }
    }
    else if ( kind == 7 )
    {
      ab = alpha;
      zemu = 2.0 / ( ab + 1.0 );

      for ( i = 0; i < m; i++ )
      {
        aj[i] = 0.0;
      }

      for ( i = 1; i <= m; i++ )
      {
        abi = i + ab * ( i % 2 );
        abj = 2 * i + ab;
        bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
      }
    }
    else if ( kind == 8 )
    {
      ab = alpha + beta;
      zemu = std::tgamma ( alpha + 1.0 ) * std::tgamma ( - ( ab + 1.0 ) ) 
        / std::tgamma ( - beta );
      apone = alpha + 1.0;
      aba = ab * apone;
      aj[0] = - apone / ( ab + 2.0 );
      bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
      for ( i = 2; i <= m; i++ )
      {
        abti = ab + 2.0 * i;
        aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
        aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
      }

      for ( i = 2; i <= m - 1; i++ )
      {
        abti = ab + 2.0 * i;
        bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) 
          / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
      }
      bj[m-1] = 0.0;
      for ( i = 0; i < m; i++ )
      {
        bj[i] =  sqrt ( bj[i] );
      }
    }

    std::cout << "leave class matrix\n" << std::endl;//DEBUG

    return zemu;
  }

  void imtqlx (int n, double d[], double e[], double z[])
  {
    //  Purpose:
    //    IMTQLX diagonalizes a symmetric tridiagonal matrix.
    //
    //  Discussion:
    //    This routine is a slightly modified version of the EISPACK routine to 
    //    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
    //
    //    The authors thank the authors of EISPACK for permission to use this
    //    routine. 
    //
    //    It has been modified to produce the product Q' * Z, where Z is an input 
    //    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
    //    The changes consist (essentialy) of applying the orthogonal transformations
    //    directly to Z as they are generated.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    08 January 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //    Roger Martin, James Wilkinson,
    //    The Implicit QL Algorithm,
    //    Numerische Mathematik,
    //    Volume 12, Number 5, December 1968, pages 377-383.
    //
    //  Parameters:
    //    Input, int N, the order of the matrix.
    //
    //    Input/output, double D(N), the diagonal entries of the matrix.
    //    On output, the information in D has been overwritten.
    //
    //    Input/output, double E(N), the subdiagonal entries of the 
    //    matrix, in entries E(1) through E(N-1).  On output, the information in
    //    E has been overwritten.
    //
    //    Input/output, double Z(N).  On input, a vector.  On output,
    //    the value of Q' * Z, where Q is the matrix that diagonalizes the
    //    input symmetric tridiagonal matrix.

    std::cout << "imtqlx\n" << std::endl;//DEBUG

    double b;
    double c;
    double f;
    double g;
    int i;
    int ii;
    int itn = 30;
    int j;
    int k;
    int l;
    int m;
    int mml;
    double p;
    double prec;
    double r;
    double s;

    prec = r8_epsilon ( );

    if ( n == 1 )
    {
      return;
    }

    e[n-1] = 0.0;

    for ( l = 1; l <= n; l++ )
    {
      j = 0;
      for ( ; ; )
      {
        for ( m = l; m <= n; m++ )
        {
          if ( m == n )
          {
            break;
          }

          if ( std::fabs ( e[m-1] ) <= prec * ( std::fabs ( d[m-1] ) + std::fabs ( d[m] ) ) )
          {
            break;
          }
        }
        p = d[l-1];
        if ( m == l )
        {
          break;
        }
        if ( itn <= j )
          Rcpp::stop("IMTQLX: Fatal error.\n\tIteration limit exceeded.");
        j = j + 1;
        g = ( d[l] - p ) / ( 2.0 * e[l-1] );
        r =  sqrt ( g * g + 1.0 );
        g = d[m-1] - p + e[l-1] / ( g + std::fabs ( r ) * r8_sign ( g ) );
        s = 1.0;
        c = 1.0;
        p = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          i = m - ii;
          f = s * e[i-1];
          b = c * e[i-1];

          if ( std::fabs ( g ) <= std::fabs ( f ) )
          {
            c = g / f;
            r =  sqrt ( c * c + 1.0 );
            e[i] = f * r;
            s = 1.0 / r;
            c = c * s;
          }
          else
          {
            s = f / g;
            r =  sqrt ( s * s + 1.0 );
            e[i] = g * r;
            c = 1.0 / r;
            s = s * c;
          }
          g = d[i] - p;
          r = ( d[i-1] - g ) * s + 2.0 * c * b;
          p = s * r;
          d[i] = g + p;
          g = c * r - b;
          f = z[i];
          z[i] = s * z[i-1] + c * f;
          z[i-1] = c * z[i-1] - s * f;
        }
        d[l-1] = d[l-1] - p;
        e[l-1] = g;
        e[m-1] = 0.0;
      }
    }
    //
    //  Sorting.
    //
    for ( ii = 2; ii <= m; ii++ )
    {
      i = ii - 1;
      k = i;
      p = d[i-1];

      for ( j = ii; j <= n; j++ )
      {
        if ( d[j-1] < p )
        {
          k = j;
          p = d[j-1];
        }
      }

      if ( k != i )
      {
        d[k-1] = d[i-1];
        d[i-1] = p;
        p = z[i-1];
        z[i-1] = z[k-1];
        z[k-1] = p;
      }
    }

    std::cout << "leave imtqlx\n" << std::endl;//DEBUG

    return;

  }

  void parchk (int kind, int m, double alpha, double beta)
  {
    //  Purpose:
    //    PARCHK checks parameters ALPHA and BETA for classical weight functions. 
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    07 January 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //  Parameters:
    //    Input, int KIND, the rule.
    //    1, Legendre,             (a,b)       1.0
    //    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
    //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    //    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    //    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    //    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
    //
    //    Input, int M, the order of the highest moment to
    //    be calculated.  This value is only needed when KIND = 8.
    //
    //    Input, double ALPHA, BETA, the parameters, if required
    //    by the value of KIND.

    std::cout << "parchk\n" << std::endl;//DEBUG

    double tmp;

    if ( kind <= 0 )
      Rcpp::stop("PARCHK: Fatal error\n\tKIND <= 0");

    //  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
    if ( 3 <= kind && alpha <= -1.0 )
      Rcpp::stop("PARCHK: Fatal error\n\t3 <= KIND and ALPHA <= -1");

    //  Check BETA for Jacobi.
    if ( kind == 4 && beta <= -1.0 )
      Rcpp::stop("PARCHK: Fatal error\n\tKIND == 4 and BETA <= -1.0");

    //  Check ALPHA and BETA for rational.
    if ( kind == 8 )
    {
      tmp = alpha + beta + m + 1.0;
      if ( 0.0 <= tmp || tmp <= beta )
        Rcpp::stop("PARCHK: Fatal error\n\tKIND == 8 but condition on ALPHA and BETA fails");
    }

    std::cout << "leave parchk\n" << std::endl;//DEBUG

    return;
  }

  double r8_epsilon (void)
  {
    //  Purpose:
    //    R8_EPSILON returns the R8 roundoff unit.
    //
    //  Discussion:
    //    The roundoff unit is a number R which is a power of 2 with the
    //    property that, to the precision of the computer's arithmetic,
    //      1 < 1 + R
    //    but
    //      1 = ( 1 + R / 2 )
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //    01 September 2012
    //
    //  Author:
    //    John Burkardt
    //
    //  Parameters:
    //    Output, double R8_EPSILON, the R8 round-off unit.

    const double value = 2.220446049250313E-016;
    return value;
  }

  double r8_sign (double x)
  {
    //  Purpose:
    //    R8_SIGN returns the sign of an R8.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    18 October 2004
    //
    //  Author:
    //    John Burkardt
    //
    //  Parameters:
    //    Input, double X, the number whose sign is desired.
    //    Output, double R8_SIGN, the sign of X.

    double value;

    if ( x < 0.0 )
    {
      value = -1.0;
    } 
    else
    {
      value = 1.0;
    }
    return value;
  }

  void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
      double swts[], double st[], int kind, double alpha, double beta, double a, 
      double b )
  {
    //  Purpose:
    //    SCQF scales a quadrature formula to a nonstandard interval.
    //
    //  Discussion:
    //    The arrays WTS and SWTS may coincide.
    //
    //    The arrays T and ST may coincide.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    16 February 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //  Parameters:
    //    Input, int NT, the number of knots.
    //    Input, double T[NT], the original knots.
    //    Input, int MLT[NT], the multiplicity of the knots.
    //    Input, double WTS[NWTS], the weights.
    //    Input, int NWTS, the number of weights.
    //    Input, int NDX[NT], used to index the array WTS. For more details see the comments in CAWIQ.
    //    Output, double SWTS[NWTS], the scaled weights.
    //    Output, double ST[NT], the scaled knots.
    //    Input, int KIND, the rule.
    //    1, Legendre,             (a,b)       1.0
    //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
    //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
    //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
    //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
    //    Input, double ALPHA, the value of Alpha, if needed.
    //    Input, double BETA, the value of Beta, if needed.
    //    Input, double A, B, the interval endpoints.

    std::cout << "scqf\n" << std::endl;//DEBUG

    double al;
    double be;
    int i;
    int k;
    int l;
    double p;
    double shft;
    double slp;
    double temp;
    double tmp;

    temp = r8_epsilon ( );

    parchk ( kind, 1, alpha, beta );

    if ( kind == 1 )
    {
      al = 0.0;
      be = 0.0;
      if ( std::fabs ( b - a ) <= temp )
        Rcpp::stop("SCQF: Fatal error\n\t|B - A| too small");
      shft = ( a + b ) / 2.0;
      slp = ( b - a ) / 2.0;
    }
    else if ( kind == 2 )
    {
      al = -0.5;
      be = -0.5;
      if ( std::fabs ( b - a ) <= temp )
        Rcpp::stop("SCQF: Fatal error\n\t|B - A| too small");
      shft = ( a + b ) / 2.0;
      slp = ( b - a ) / 2.0;
    }
    else if ( kind == 3 )
    {
      al = alpha;
      be = alpha;
      if ( std::fabs ( b - a ) <= temp )
        Rcpp::stop("SCQF: Fatal error\n\t|B - A| too small");
      shft = ( a + b ) / 2.0;
      slp = ( b - a ) / 2.0;
    }
    else if ( kind == 4 )
    {
      al = alpha;
      be = beta;

      if ( std::fabs ( b - a ) <= temp )
        Rcpp::stop("SCQF: Fatal error\n\t|B - A| too small");
      shft = ( a + b ) / 2.0;
      slp = ( b - a ) / 2.0;
    }
    else if ( kind == 5 )
    {
      if ( b <= 0.0 )
        Rcpp::stop("SCQF: Fatal error\n\tB <= 0");
      shft = a;
      slp = 1.0 / b;
      al = alpha;
      be = 0.0;
    }
    else if ( kind == 6 )
    {
      if ( b <= 0.0 )
        Rcpp::stop("SCQF: Fatal error\n\tB <= 0");
      shft = a;
      slp = 1.0 / sqrt ( b );
      al = alpha;
      be = 0.0;
    }
    else if ( kind == 7 )
    {
      al = alpha;
      be = 0.0;
      if ( std::fabs ( b - a ) <= temp )
        Rcpp::stop("SCQF: Fatal error\n\t|B - A| too small");
      shft = ( a + b ) / 2.0;
      slp = ( b - a ) / 2.0;
    }
    else if ( kind == 8 )
    {
      if ( a + b <= 0.0 )
        Rcpp::stop("SCQF: Fatal error\n\tA + B <= 0");
      shft = a;
      slp = a + b;
      al = alpha;
      be = beta;
    }
    else if ( kind == 9 )
    {
      al = 0.5;
      be = 0.5;
      if ( std::fabs ( b - a ) <= temp )
        Rcpp::stop("SCQF: Fatal error\n\t|B - A| too small");
      shft = ( a + b ) / 2.0;
      slp = ( b - a ) / 2.0;
    }

    p = std::pow ( slp, al + be + 1.0 );

    for ( k = 0; k < nt; k++ )
    {
      st[k] = shft + slp * t[k];
      l = std::abs ( ndx[k] );

      if ( l != 0 )
      {
        tmp = p;
        for ( i = l - 1; i <= l - 1 + mlt[k] - 1; i++ )
        {
          swts[i] = wts[i] * tmp;
          tmp = tmp * slp;
        }
      }
    }

    std::cout << "leave scqf\n" << std::endl;//DEBUG
    return;
  }

  void sgqf (int nt, double aj[], double bj[], double zemu, double t[], double wts[])
  {
    //  Purpose:
    //    SGQF computes knots and weights of a Gauss Quadrature formula.
    //
    //  Discussion:
    //    This routine computes all the knots and weights of a Gauss quadrature
    //    formula with simple knots from the Jacobi matrix and the zero-th
    //    moment of the weight function, using the Golub-Welsch technique.
    //
    //  Licensing:
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //    08 January 2010
    //
    //  Author:
    //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //    Sylvan Elhay, Jaroslav Kautsky,
    //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    //    Interpolatory Quadrature,
    //    ACM Transactions on Mathematical Software,
    //    Volume 13, Number 4, December 1987, pages 399-415.
    //
    //  Parameters:
    //    Input, int NT, the number of knots.
    //    Input, double AJ[NT], the diagonal of the Jacobi matrix.
    //    Input/output, double BJ[NT], the subdiagonal of the Jacobi matrix, in entries 1 through NT-1. On output, BJ has been overwritten.
    //    Input, double ZEMU, the zero-th moment of the weight function.
    //    Output, double T[NT], the knots.
    //    Output, double WTS[NT], the weights.

    std::cout << "sgqf\n" << std::endl;//DEBUG

    int i;

    //  Exit if the zero-th moment is not positive.
    if ( zemu <= 0.0 )
      Rcpp::stop("SGQF: Fatal error");

    //  Set up vectors for IMTQLX.
    for ( i = 0; i < nt; i++ )
    {
      t[i] = aj[i];
    }
    wts[0] = sqrt ( zemu );
    for ( i = 1; i < nt; i++ )
    {
      wts[i] = 0.0;
    }

    //  Diagonalize the Jacobi matrix.
    imtqlx ( nt, t, bj, wts );

    for ( i = 0; i < nt; i++ )
    {
      wts[i] = wts[i] * wts[i];
    }

    std::cout << "leave sgqf\n" << std::endl;//DEBUG

    return;
  }
};

//--------------------------------------- wrappers around specific rules

struct LaguerreRule : public GaussianQuadratureRules
{
  // WHY
  //   Quadrature rule for integrating functions of the form g(x) = x^{\alpha-1} e^{-\beta x} over [0, inf].
  //   This is equivalent to marginalizing with regards to a gamma distribution.
  //   The normalizing constant is \gamma(\alpha)/\beta^\alpha
  // HOW
  //   order : order of rule
  //   alpha : real; alpha > 0 
  //   beta  : real; beta > 0
  // WHAT
  //   abcissa, weights stored as members

  uword dim, order;
  mat abscissa; 
  vec weights;

  LaguerreRule (uword order, double alpha = 1.0, double beta = 1.0)
    : order(order)
    , dim(1)
    , weights(order)
    , abscissa(order, 1)
  {
    cgqf (order, 5, alpha - 1.0, 0.0, 0.0, beta, abscissa.memptr(), weights.memptr());
    double constant = std::pow(beta, alpha) / std::tgamma(alpha);
    weights *= constant;
  }
};

struct HermiteRule : public GaussianQuadratureRules
{
  // WHY
  //   Quadrature rule for integrating functions of the form g(x) = \sqrt{1/(2\pi)} e^{-x^2} f(\mu + x/\sqrt{tau}) over [-inf, inf]
  //   This is equivalent to marginalizing with regards to a univariate Gaussian distribution.
  // HOW
  //   order : order of rule
  //   mu    : real; -inf < center < inf
  //   tau   : real; tau > 0
  // WHAT
  //   abcissa, weights stored as members

  uword dim, order;
  mat abscissa; 
  vec weights;

  HermiteRule (uword order, double mu = 0.0, double tau = 1.0)
    : order(order)
    , dim(1)
    , weights(order)
    , abscissa(order, 1)
  {
    cgqf (order, 6, 0.0, 0.0, mu, tau, abscissa.memptr(), weights.memptr());
    double constant = std::sqrt(tau / arma::datum::pi);
    weights *= constant;
  }
};

struct HomeierRule : public GaussianQuadratureRules
{
  // WHY
  //   Quadrature rule for integrating functions of the form g(x) = 2\gamma/\pi \frac{1}{x^{0.5}(x+\gamma^2)} f(x^{0.5}) over [0, \inf].
  //   This is equivalent to marginalizing with regards to a half-Cauchy distribution with scale \gamma.
  //
  //   The core idea is to use a Mobius transformation, x = \gamma^2 (1-u)/(1+u), and then integrate the function
  //
  //     1/\pi (1 + u)^{-0.5} (1 - u)^{-0.5} f(gamma/\sqrt{(1 - u)}\sqrt{(1 + u)} )
  //
  //   over [-1, 1] using a using a Gauss-Chebyshev rule. Reference is:
  //
  //     Homeier HHH, Steinborn EO. 1996. Computer Physics Communications 99: 77-80
  //   
  // HOW
  //   order : order of rule
  //   gamma : real; gamma > 0
  // WHAT
  //   abcissa, weights stored as members

  uword dim, order;
  mat abscissa; 
  vec weights;

  HomeierRule (uword order, double gamma = 1.0)
    : order(order)
    , dim(1)
    , weights(order, arma::fill::zeros)
    , abscissa(order, 1, arma::fill::zeros)
  {
    if (order < 1 || gamma <= 0.0)
      Rcpp::stop("HomeierRule: invalid parameters");
    cgqf (order, 3, -0.5, 0.0, -1, 1, abscissa.memptr(), weights.memptr());
    weights /= arma::datum::pi;
    abscissa = gamma * arma::sqrt((1.+abscissa)/(1.-abscissa));
  }
};

struct CauchyRule
{
  // weight: all equal to 1/order
  // abscissa: x = cos((2*i - 1)/(2*order)*pi), gamma*sqrt((1-x)/(1+x))

};

struct JacobiRule : public GaussianQuadratureRules
{
  // TODO
};

struct LegendreRule : public GaussianQuadratureRules
{
  // TODO
};

struct ChebyshevRule : public GaussianQuadratureRules
{
  // TODO
};

struct ExponentialRule : public GaussianQuadratureRules
{
  // TODO
};

struct RationalRule : public GaussianQuadratureRules
{
  // TODO
};

//---------------------- tests

//[[Rcpp::export("inlassle_test_HomeierRule")]]
Rcpp::List test_HomeierRule (arma::uword order, double gamma)
{
  HomeierRule rule (order, gamma);
  return Rcpp::List::create (
      Rcpp::_["abscissa"] = rule.abscissa,
      Rcpp::_["weights"] = rule.weights);
}
