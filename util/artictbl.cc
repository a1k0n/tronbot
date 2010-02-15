// build articulation table for O(1) updates of connected-component algo in the
// best cases

#include <stdio.h>

bool is_open(int n, int min, int max)
{
  for(int i=min;i<=max;i++)
    if(n&(1<<(i&7))) return false;
  return true;
}

int potential_articulation(int n)
{
  int degree = 0;
  for(int i=1;i<7;i+=2) {
    if(n&(1<<i)) continue;
    degree++;
    for(int j=i+2;j<=7;j+=2) {
      if(n&(1<<j)) continue;
      if(!is_open(n, i+1, j-1) && !is_open(n, j+1, i+7)) {
//        printf("%d is closed between %d and %d\n", n, i, j);
        return 1;
      }
    }
  }
  return 0;
}

int main()
{
  printf("char _potential_articulation[256] = {\n  ");
  for(int n=0;n<256;n++) {
    printf("%s%d", n==0 ? "" : (n&15) == 0 ? ",\n  " : ",",
           potential_articulation(n));
  }
}

