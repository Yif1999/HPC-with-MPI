#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MAX 10
typedef double dtype;
dtype num_genarator (int i,char *s);

void main (int argc,char *argv[]){
  char *mode=argv[1];
  int r=atoi(argv[2]),c=atoi(argv[3]);
  char *filename=argv[4];
  int count=0;
  int i,j;
  FILE  *fileptr; 
  fileptr=fopen(filename,"w");

  fwrite(&r, sizeof(int), 1, fileptr);
  // fwrite(&c, sizeof(int), 1, fileptr);
  for (i=0;i<r;i++){
    for (j=0;j<c;j++){
      count++;
      dtype num=num_genarator(count,mode);
      fwrite(&num, sizeof(dtype), 1, fileptr);
    }
  }
  fclose(fileptr);

  return;
}

dtype num_genarator (int i,char *mode){
  if (strcmp(mode,"random")==0){
    srand((unsigned)time(NULL));
    return rand()%10+1;
    printf("Random mat generated successfully!\n");
  }
  else if (strcmp(mode,"order")==0){
    return i;
    printf("Order mat generated successfully!\n");
  }
  else{
    printf("[ERROR] Invalid Input!\n");
  }
  return 0;
}
