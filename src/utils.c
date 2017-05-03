#include "utils.h"

int cmp_int(void *a, void *b)
{
  if(*(int*)a > *(int*)b)
  {
    return 1;
  }else if(*(int*)a == *(int*)b)
  {
    return 0;
  }else{
    return -1;
  }
}

int cmp_dou(void *a, void *b)
{
  if(*(double*)a > *(double*)b)
  {
    return 1;
  }else if(*(double*)a == *(double*)b)
  {
    return 0;
  }else{
    return -1;
  }
}

int strBin2Dec(char *str, int *dec){

  /* Converts the string in "str" that suppose to have a continuue
   * sequence of "0" and "1" to its decimal representation.
   */

  int i;
  *dec=0; 
  for(i=strlen(str)-1;i>=0;i--){
    if(str[i]=='0' || str[i]=='1'){
      if(str[i]=='1')
        *dec+=(int)pow(2,strlen(str)-1-i);
    }else{
      return 1;
    }
  }
  return 0;
}
