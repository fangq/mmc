#include <stdio.h>
#include <smmintrin.h>

void main(){
	__m128 O,T,S;
        float len;
	const float tmax=1e10f;
	int idx;
	const char maskmap[16]={4,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3};

	T = _mm_set_ps(0.01f,0.1f,0.002f,0.1f);
	S = _mm_set_ps(1.0, 0.0, -1.0f, 1.0);

	O = _mm_cmpge_ps(S,_mm_set1_ps(0.f));
	T = _mm_add_ps(_mm_andnot_ps(O,_mm_set1_ps(tmax)),_mm_and_ps(O,T));
	S = _mm_movehl_ps(T, T);
	O = _mm_min_ps(T, S);
	S = _mm_shuffle_ps(O, O, _MM_SHUFFLE(1,1,1,1));
	O = _mm_min_ss(O, S);
	
	_mm_store_ss(&len,O);
	idx=_mm_movemask_ps(_mm_cmpeq_ps(T,_mm_set1_ps(len)));
	printf("%f %d %d\n",len,idx,maskmap[idx]);
}
