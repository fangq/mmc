#include <set>
#include <list>
#include <algorithm>
#include <string.h>
#include "simpmesh.h"
#include "highordermesh.h"
#include "mcx_utils.h"
#include <unordered_map>

#define TETEDGE 6

const int edgepair[TETEDGE][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

#ifdef __cplusplus
extern "C"
#endif

struct Key
{
  public:
  int first;
  int second;
  bool operator==(const Key &other) const
  { return (first == other.first
            && second == other.second);
  }
};

namespace std {

  template <>
  struct hash<Key>{
    std::size_t operator()(const Key& k) const {
      using std::size_t;
      using std::hash;

      if(sizeof(std::size_t)==8)
          return ( ((std::size_t)(hash<int>()(k.first) ))<< 32
                | (hash<int>()(k.second) ) );
      else
          return ((hash<int>()(k.first)
               ^ (hash<int>()(k.second) << 1)) >> 1);
    }
  };

}

void mesh_10nodetet(tetmesh * mesh,mcconfig *cfg){
	int pos,n1,n2,oldnn=mesh->nn;
	int *ee;
	int Count=mesh->nn; // if from 1 to nn -> Count=nn+1

	std::unordered_map<Key,int> edgemap;
	std::unordered_map<Key,int>::const_iterator imap;
  	std::list<std::pair<int,int> > edgelist;
  	std::list<std::pair<int,int> >::iterator it;

	std::pair<int,int> edge;
	it = edgelist.begin();

	if(mesh->elem2==NULL)
	    mesh->elem2=(int *)calloc(sizeof(int)*TETEDGE,mesh->ne);

	mesh->nvol=(float *)realloc((void*)mesh->nvol,sizeof(float)*(mesh->nn+mesh->ne*TETEDGE));
	memset(mesh->nvol+mesh->nn,0,sizeof(float)*mesh->ne*TETEDGE);

	for (int eid=0;eid<mesh->ne;eid++)
	    for(int ed=0; ed<TETEDGE; ed++){
	        ee=(int*)(&mesh->elem[eid]);
		n1=MIN(ee[edgepair[ed][0]],ee[edgepair[ed][1]])-1;
		n2=MAX(ee[edgepair[ed][0]],ee[edgepair[ed][1]])-1;
		edge=std::make_pair(n1,n2);

		imap = edgemap.find ({n1,n2});//find
		if ( imap == edgemap.end() ){ // no such edge was found
			edgemap.insert (std::make_pair<Key,int>({n1,n2},std::add_rvalue_reference<int>::type(Count+1)));//insert
			edgelist.insert(it,edge);
			Count+=1;// what about count?
                        mesh->elem2[eid*TETEDGE+ed]=Count;
		}else{ // an existing edge was found in the hash table
                        mesh->elem2[eid*TETEDGE+ed]=imap->second;
                }
		mesh->nvol[mesh->elem2[eid*TETEDGE+ed]-1]+=mesh->evol[eid];
	}
	pos=0;

	mesh->nn+=edgelist.size();
	mesh->nvol=(float *)realloc((void*)mesh->nvol,sizeof(float)*(mesh->nn));
	for(int i=0;i<oldnn;i++)
		mesh->nvol[i]=-mesh->nvol[i]*(0.05f*4.f); // 4 undo the 1st order normalization, 0.05 is the normalization for the 2nd basis on the corner nodes
	for(int i=oldnn;i<mesh->nn;i++)
		mesh->nvol[i]*=0.2f;
	mesh->node=(float3*)realloc((void*)mesh->node, sizeof(float3)*(mesh->nn));
        mesh->weight=(double *)realloc((void*)mesh->weight,sizeof(double)*mesh->nn*cfg->maxgate);
        memset(mesh->weight,0,sizeof(double)*mesh->nn*cfg->maxgate); // if mesh->weight is filled, need to allocate a new buffer, and copy the old buffer gate by gate

	for (it=edgelist.begin(); it!=edgelist.end(); it++){
	    for(int i=0;i<3;i++){
	        ((float*)(&mesh->node[oldnn+pos]))[i]= 
		   (((float*)(&mesh->node[(*it).first]))[i]+((float*)(&mesh->node[(*it).second]))[i])*0.5f;
	    }
	    pos++;
	}
}
