using namespace std;

#include <iostream>
#include "randomGen.h"

randomGen sR;

int main (void)
{
  double sum, tmp;
  int cnt1= 0, cnt2=0;
  
  sum= 0.0;
  for (int i= 0; i < 10000; i++) {
    tmp= sR.n();
    sum+= tmp;
    if (tmp < 0.4) cnt1++;
    else cnt2++;
    cerr << tmp << " " << sum << endl;
  }
  cerr << cnt1 << " " << cnt2 << endl;
  return 0;
}
