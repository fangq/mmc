#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef MMC_LOGISTIC
  #include "logistic_rand.c"
#elif defined MMC_SFMT
  #include "sfmt_rand.c"
#else
  #include "posix_randr.c"
#endif

void usage(char *exename){
	printf("usage:\n\t%s <count|100000> <seed|19650218>\n",exename);
}
void genrand(int threadid, int id, RandType ran[],RandType ran0[]){
	rand_need_more(ran,ran0);
	printf("%d\t%d\t%e\t%e\t%e\n",threadid,id,rand_next_scatlen(ran),rand_next_aangle(ran),rand_next_zangle(ran));
}
int main(int argc, char**argv){
	unsigned int seed=19650218, count=100000,threadid=0,i;
        RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
        RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));

	if(argc>=2)
		count=atoi(argv[1]);
        if(argc>=3) 
                seed=atoi(argv[2]);

	if(count==0 || seed==0 || argc==1){ /*when atoi returned error*/
		usage(argv[0]);
		exit(0);
	}

	if(seed<0) seed=time(NULL);

#pragma omp parallel private(ran0,ran1,threadid)
{
#ifdef _OPENMP
	threadid=omp_get_thread_num();	
#endif
	rng_init(ran0,ran1,(unsigned int *)&seed,threadid);

#pragma omp for
	for(i=0;i<count;i++){
		genrand(threadid,i,ran0,ran1);
	}
}
	return 0;
}
