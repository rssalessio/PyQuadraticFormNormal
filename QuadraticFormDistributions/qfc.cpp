#include <cmath>
#include <malloc.h>
#include <setjmp.h>

using std::atan;
using std::exp;
using std::sin;
using std::log;
using std::fabs;
using std::sqrt;
using std::floor;


#define PI 3.14159265358979323846264338328
#define log28 .0866  /*  log(2.0) / 8.0  */
 

/* to avoid underflows  */
static double exp1(double x) {
   return x < -50.0 ? 0.0 : std::exp(x);
}

static double square(double x)  {
   return x*x;
}

static double cube(double x)  {
   return x*x*x;
}

/* if (first) log(1 + x) ; else  log(1 + x) - x */
static double  log1(double x, bool first) {
   if (std::fabs(x) > 0.1) {
      return (first ? std::log(1.0 + x) : (std::log(1.0 + x) - x));
   } else {
      double s, s1, term, y, k;
      y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
      s = (first ? 2.0 : - x) * y;
      y = square(y);
      for (s1 = s + term / k; s1 != s; s1 = s + term / k)
      { k = k + 2.0; term = term * y; s = s1; }
      return s;
   }
}


struct DaviesMethodParameters {
   DaviesMethodParameters(
         const double *lambda, const double *nc, const ssize_t *df,
         const ssize_t limit, const ssize_t r, double x, double sigma): 
         lambda(lambda), nc(nc), df(df), limit(limit), counter(0), r(r), ndtsrt(true), fail(false),
         sigsq(square(sigma)), lmax(0.0), lmin(0.0), mean(0), c(x), intl(0), ersm(0), sigma(sigma) {
      th = new ssize_t[r];
   }
   ~DaviesMethodParameters() {
      delete []th;
   }

   /*  count number of calls to errbd, truncation, cfe */
   bool step_counter() {
      return ++counter  > limit;
   }

   const double *lambda, *nc;
   const ssize_t *df;
   ssize_t counter, *th;
   const ssize_t limit, r;
   bool ndtsrt, fail;
   jmp_buf env;
   double sigsq, lmax, lmin, mean, c, intl, ersm, sigma;
};


/* find order of absolute values of lb */
static void order(DaviesMethodParameters &parameters) {
   int j, k; double lj;
   for ( j = 0; j < parameters.r; j++ ) {
      lj = std::fabs(parameters.lambda[j]);
      for (k = j - 1; k >= 0; k--) {
         if ( lj > std::fabs(parameters.lambda[parameters.th[k]]) )
            parameters.th[k + 1] = parameters.th[k];
         else
            goto l1;
      }
      k = -1;
   l1 :
      parameters.th[k + 1] = j;
   }
   parameters.ndtsrt = false;
}


/*  find bound on tail probability using mgf, cutoff
      point returned to *cx */
static double errbd(double u, double* cx, DaviesMethodParameters &parameters) {
   double sum1, xconst;
   
   if (parameters.step_counter())
      longjmp(parameters.env, 1);

   xconst = u * parameters.sigsq; 
   sum1 = u * xconst;
   u = 2.0 * u;
   
   for (int j=parameters.r-1; j>=0; j--) {
      ssize_t nj = parameters.df[j];
      double lj = parameters.lambda[j];
      double ncj = parameters.nc[j];
      double x = u * lj;
      double y = 1.0 - x;
      xconst = xconst + lj * (ncj / y + nj) / y;
      sum1 = sum1 + ncj * square(x / y)
         + nj * (square(x) / y + log1(-x, false ));
   }
   *cx = xconst;
   return exp1(-0.5 * sum1);
}

/*  find ctff so that p(qf > ctff) < accx  if (upn > 0,
       p(qf < ctff) < accx otherwise */
static double  ctff(double accx, double* upn, DaviesMethodParameters &parameters) {
   double u1, u2, u, rb, xconst, c1, c2;
   u2 = *upn;   u1 = 0.0;  c1 = parameters.mean;
   rb = 2.0 * ((u2 > 0.0) ? parameters.lmax : parameters.lmin);
   for (u = u2 / (1.0 + u2 * rb); errbd(u, &c2, parameters) > accx; 
      u = u2 / (1.0 + u2 * rb)) {
      u1 = u2;
      c1 = c2;
      u2 = 2.0 * u2;
   }
   for (u = (c1 - parameters.mean) / (c2 - parameters.mean); u < 0.9;
      u = (c1 - parameters.mean) / (c2 - parameters.mean)) {
      u = (u1 + u2) / 2.0;
      if (errbd(u / (1.0 + u * rb), &xconst, parameters) > accx) { 
         u1 = u;
         c1 = xconst;
      }
      else {
         u2 = u;
         c2 = xconst;
      }
   }
   *upn = u2;
   return c2;
}

/* bound integration error due to truncation at u */
static double truncation(double u, double tausq, DaviesMethodParameters &parameters) {
   double sum1, sum2, prod1, prod2, prod3, lj, ncj,
      x, y, err1, err2;
   int j, nj, s;
   if (parameters.step_counter())
      longjmp(parameters.env, 1);
   sum1  = 0.0; prod2 = 0.0;  prod3 = 0.0;  s = 0;
   sum2 = (parameters.sigsq + tausq) * square(u); prod1 = 2.0 * sum2;
   u = 2.0 * u;
   for (j=0; j<parameters.r; j++ ) {
      lj = parameters.lambda[j];  ncj = parameters.nc[j]; nj = parameters.df[j];
      x = square(u * lj);
      sum1 = sum1 + ncj * x / (1.0 + x);
      if (x > 1.0) {
         prod2 = prod2 + nj * log(x);
         prod3 = prod3 + nj * log1(x, true );
         s = s + nj;
      }
      else  prod1 = prod1 + nj * log1(x, true );
   }
   sum1 = 0.5 * sum1;
   prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
   x = exp1(-sum1 - 0.25 * prod2) / PI;
   y = exp1(-sum1 - 0.25 * prod3) / PI;
   err1 =  ( s  ==  0 )  ? 1.0 : x * 2.0 / s;
   err2 =  ( prod3 > 1.0 )  ? 2.5 * y : 1.0;
   if (err2 < err1)
      err1 = err2;
   x = 0.5 * sum2;
   err2 =  ( x  <=  y )  ? 1.0  : y / x;
   return  ( err1 < err2 )  ? err1  :  err2;
}

/*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
static void findu(double* utx, double accx, DaviesMethodParameters &parameters) {
   double u, ut; int i;
   static double divis[]={2.0,1.4,1.2,1.1};
   ut = *utx; u = ut / 4.0;
   if ( truncation(u, 0.0, parameters) > accx ) {
      for ( u = ut; truncation(u, 0.0, parameters) > accx; u = ut)
         ut = ut * 4.0;
   }
   else {
      ut = u;
      for( u = u / 4.0; truncation(u, 0.0, parameters) <=  accx; u = u / 4.0 )
         ut = u;
   }
   for ( i=0;i<4;i++) {
      u = ut/divis[i];
      if ( truncation(u, 0.0, parameters)  <=  accx )
         ut = u;
   }
   *utx = ut;
}


/*  carry out integration with nterm terms, at stepsize
   interv.  if (! mainx) multiply integrand by
      1.0-exp(-0.5*tausq*u^2) */
static void integrate(int nterm, double interv, double tausq, bool mainx, DaviesMethodParameters &parameters) {
   double inpi, u, sum1, sum2, sum3, x, y, z;
   int k, j, nj;
   inpi = interv / PI;
   for ( k = nterm; k>=0; k--) {
      u = (k + 0.5) * interv;
      sum1 = - 2.0 * u * parameters.c;
      sum2 = std::fabs(sum1);
      sum3 = - 0.5 * parameters.sigsq * square(u);
      for ( j = parameters.r-1; j>=0; j--) {
         nj = parameters.df[j];
         x = 2.0 * parameters.lambda[j] * u;
         y = square(x);
         sum3 = sum3 - 0.25 * nj * log1(y, true );
         y = parameters.nc[j] * x / (1.0 + y);
         z = nj * std::atan(x) + y;
         sum1 = sum1 + z;
         sum2 = sum2 + std::fabs(z);
         sum3 = sum3 - 0.5 * x * y;
      }
      x = inpi * exp1(sum3) / u;
      if (!mainx)
         x = x * (1.0 - exp1(-0.5 * tausq * square(u)));
      sum1 = std::sin(0.5 * sum1) * x;
      sum2 = 0.5 * sum2 * x;
      parameters.intl = parameters.intl + sum1;
      parameters.ersm = parameters.ersm + sum2;
   }
}

/*  coef of tausq in error when convergence factor of
   exp1(-0.5*tausq*u^2) is used when df is evaluated at x */

static double cfe(double x, DaviesMethodParameters &parameters) {
   double axl, axl1, axl2, sxl, sum1, lj;
   int j, k, t;
   if (parameters.step_counter())
      longjmp(parameters.env, 1);
   if (parameters.ndtsrt)
      order(parameters);
   axl = std::fabs(x);  sxl = (x>0.0) ? 1.0 : -1.0;  sum1 = 0.0;
   for ( j = parameters.r-1; j>=0; j-- ){
      t = parameters.th[j];
      if ( parameters.lambda[t] * sxl > 0.0 ) {
         lj = std::fabs(parameters.lambda[t]);
         axl1 = axl - lj * (parameters.df[t] + parameters.nc[t]);
         axl2 = lj / log28;
         if ( axl1 > axl2 )
            axl = axl1;
         else {
            if ( axl > axl2 )  axl = axl2;
            sum1 = (axl - axl1) / lj;
            for ( k = j-1; k>=0; k--)
            sum1 = sum1 + (parameters.df[parameters.th[k]] + parameters.nc[parameters.th[k]]);
            goto  l;
         }
      }
   }
l:
   if (sum1 > 100.0) {
      parameters.fail = true;
      return 1.0;
   }
   return std::pow(2.0,(sum1 / 4.0)) / (PI * square(axl));
}


double qf(double* lb1, double* nc1, ssize_t* n1, ssize_t r1, double sigma, double c1,
   ssize_t lim1, double acc, double* trace, int* ifault)

/*  distribution function of a linear combination of non-central
   chi-squared random variables :

input:
   lb[j]            coefficient of j-th chi-squared variable
   nc[j]            non-centrality parameter
   n[j]             degrees of freedom
   j = 0, 2 ... r-1
   sigma            coefficient of standard normal variable
   c                point at which df is to be evaluated
   lim              maximum number of terms in integration
   acc              maximum error

output:
   ifault = 1       required accuracy NOT achieved
            2       round-off error possibly significant
            3       invalid parameters
            4       unable to locate integration parameters
            5       out of memory

   trace[0]         absolute sum
   trace[1]         total number of integration terms
   trace[2]         number of integrations
   trace[3]         integration interval in final integration
   trace[4]         truncation point in initial integration
   trace[5]         s.d. of initial convergence factor
   trace[6]         cycles to locate integration parameters     */

{
      int j, nj, nt, ntm;  double acc1, almx, xlim, xnt, xntm;
      double utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj, ncj;
      double qfval;
      static int rats[]={1,2,4,8};

      DaviesMethodParameters parameters = DaviesMethodParameters(lb1, nc1, n1, lim1, r1, c1, sigma);
      if (setjmp(parameters.env) != 0) {
         *ifault=4;
         goto endofproc;
      }

      for ( j = 0; j<7; j++ )
         trace[j] = 0.0;

      *ifault = 0;
      qfval = -1.0; acc1 = acc;

      xlim = (double) parameters.limit;
      if (!parameters.th) {
         *ifault=5; 
         goto  endofproc;
      } 

      /* find mean, sd, max and min of lb,
         check that parameter values are valid */
      sd = parameters.sigsq;

      for (j=0; j<parameters.r; j++ ) {
         nj = parameters.df[j]; 
         lj = parameters.lambda[j]; 
         ncj = parameters.nc[j];
         if ( nj < 0  ||  ncj < 0.0 ) {
            *ifault = 3;
            goto  endofproc;
         }
         sd  = sd  + square(lj) * (2 * nj + 4.0 * ncj);
         parameters.mean = parameters.mean + lj * (nj + ncj);
         if (parameters.lmax < lj)
            parameters.lmax = lj;
         else if (parameters.lmin > lj)
            parameters.lmin = lj;
      }

      if ( sd == 0.0  ){
         qfval = (parameters.c > 0.0) ? 1.0 : 0.0;
         goto  endofproc;
      }
      if ( parameters.lmin == 0.0 && parameters.lmax == 0.0 && parameters.sigma == 0.0 ){
         *ifault = 3;
         goto  endofproc;
      }
      sd = sqrt(sd);
      almx = (parameters.lmax < - parameters.lmin) ? - parameters.lmin : parameters.lmax;

      /* starting values for findu, ctff */
      utx = 16.0 / sd;  up = 4.5 / sd;  un = - up;
      /* truncation point with no convergence factor */
      findu(&utx, .5 * acc1, parameters);
      /* does convergence factor help */
      if (parameters.c != 0.0  && (almx > 0.07 * sd))
      {
         tausq = .25 * acc1 / cfe(parameters.c, parameters);
         if (parameters.fail)
            parameters.fail = false ;
         else if (truncation(utx, tausq, parameters) < .2 * acc1)
         {
            parameters.sigsq = parameters.sigsq + tausq;
            findu(&utx, .25 * acc1, parameters);
            trace[5] = std::sqrt(tausq);
         }
      }
      trace[4] = utx;
      acc1 = 0.5 * acc1;

      /* find RANGE of distribution, quit if outside this */
   l1:
      d1 = ctff(acc1, &up, parameters) - parameters.c;
      if (d1 < 0.0) {
         qfval = 1.0;
         goto endofproc;
      }
      d2 = parameters.c - ctff(acc1, &un, parameters);
      if (d2 < 0.0) {
         qfval = 0.0;
         goto endofproc;
      }
      /* find integration interval */
      intv = 2.0 * PI / ((d1 > d2) ? d1 : d2);
      /* calculate number of terms required for main and
         auxillary integrations */
      xnt = utx / intv;  xntm = 3.0 / std::sqrt(acc1);
      if (xnt > xntm * 1.5) {
         /* parameters for auxillary integration */
         if (xntm > xlim) {
            *ifault = 1;
            goto endofproc;
         }
         ntm = (int) std::floor(xntm+0.5);
         intv1 = utx / ntm;  x = 2.0 * PI / intv1;
         if (x <= std::fabs(parameters.c)) goto l2;
         /* calculate convergence factor */
         tausq = .33 * acc1 / (1.1 * (cfe(parameters.c - x, parameters) + cfe(parameters.c + x, parameters)));
         if (parameters.fail) goto l2;
         acc1 = .67 * acc1;
         /* auxillary integration */
         integrate(ntm, intv1, tausq, false, parameters);
         xlim = xlim - xntm; 
         parameters.sigsq = parameters.sigsq + tausq;
         trace[2] = trace[2] + 1;
         trace[1] = trace[1] + ntm + 1;
         /* find truncation point with new convergence factor */
         findu(&utx, .25 * acc1, parameters);
         acc1 = 0.75 * acc1;
         goto l1;
      }

      /* main integration */
   l2:
      trace[3] = intv;
      if (xnt > xlim) {
         *ifault = 1;
         goto endofproc;
      }
      nt = (int) std::floor(xnt+0.5);
      integrate(nt, intv, 0.0, true, parameters);
      trace[2] = trace[2] + 1;
      trace[1] = trace[1] + nt + 1;
      qfval = 0.5 - parameters.intl;
      trace[0] = parameters.ersm;

      /* test whether round-off error could be significant
         allow for radix 8 or 16 machines */
      up=parameters.ersm;
      x = up + acc / 10.0;
      for (j=0;j<4;j++) {
         if (rats[j] * x == rats[j] * up)
            *ifault = 2;
      }

   endofproc :
      trace[6] = (double)parameters.counter;
      return qfval;
}


