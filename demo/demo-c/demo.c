#include <stdio.h>
#include <stdlib.h>
#include "seuif97.h"

int main(void)
{

   double in1 = 1.07;
   double in2 = 2000;
   double out;

   out = ph(in1, in2, 2);
   out = ph(in1, in2, 19);
   out = ph(in1, in2, 20);
   out = ph(in1, in2, 21);
   out = ph(in1, in2, 22);
   out = ph(in1, in2, 30);
   out = ph(in1, in2, 31);
   printf("(in1,in2)(%.5f,%.5f) r= %.10f", in1, in2, out);
   return EXIT_SUCCESS;
}
