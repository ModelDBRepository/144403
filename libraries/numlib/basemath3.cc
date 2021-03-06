#include "basemath3.h"

double binomial(double n, double k)
{
  if (k > n/2) return binomial(n, n-k);
  
  assert(n >= 0.0);
  assert(k >= 0.0);

  double x= 0.0;
  double bi= 1;

  while (k > x)
  {
    bi*= (n-x);
    bi/= (k-x);
    x++;
  }
  return bi;
}

long double binomial(long double n, long double k)
{
  if (k > n/2) return binomial(n, n-k);

  assert(n >= 0.0L);
  assert(k >= 0.0L);
  
  long double x= 0.0L;
  long double bi= 1.0L;

  while (k > x)
  {
    bi*= (n-x);
    bi/= (k-x);
    x++;
  }
  return bi;
}

int binomial(int n, int k)
{
  assert(n >= 0);
  assert(k >= 0);

  int x= 1;
  int j= max(k, n-k);
  int l= min(k, n-k);

  for (int i= n; i > j; i--)
  {
    x*= i;
  }
  for (int i= l; i > 1; i--)
  {
    x/= i;
  }
  return x;
}

template <class Type>
inline Type ipower(Type x, int p)
{
  assert(p > 0);
  Type y= x;
  for (int i= 0; i < p-1; i++) y*= x;
  return y;
}

double p_binomial(double n, double p, double k)
{
  if (fabs(p-1.0) < 1.0e-9) {
    if (fabs(k-n) < 1.0e-9) return 1.0;
    else return 0.0;
  }
  
  if (p < 1.0e-9) {
    if (k < 1.0e-9) return 1.0;
    else return 0.0;
  }

  double val= binomial(n,k);
  val*= exp(k*log(p));
  val*= exp((n-k)*log(1.0-p));
      
  return val;
}

long double p_binomial(long double n, long double p, long double k)
{
  static long double val, expo;
  
  if (fabsl(p-1.0L) < 1.0e-9L) {
    if (fabsl(k-n) < 1.0e-9L) return 1.0L;
    else return 0.0L;
  }
  
  if (p < 2.0e-9L) {
    if (k < 1.0e-9L) return 1.0L;
    else return 0.0L;
  }

  val= binomial(n,k);
  if (isinf(val)) {
    cerr << "binomial(" << n << ", " << k << ") is inf" << endl;
  }
  expo= k*logl(p);
  val*= expl(expo);
  expo= (n-k)*logl(1.0L-p);
  val*= expl(expo);

  return val;
}




