#include <set>
#include <list>
#include <algorithm>
#include <string.h>
#include "simpmesh.h"
#include "highordermesh.h"
#include "mcx_utils.h"

#define TETEDGE 6

const int edgepair[TETEDGE][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

#ifdef __cplusplus
extern "C"
#endif
void mesh_10nodetet(tetmesh * mesh,mcconfig *cfg){
	int pos,n1,n2,oldnn=mesh->nn;
	unsigned int setlen;
	int *ee;

  	std::set<std::pair<int,int> > edgeset;
  	std::set<std::pair<int,int> >::iterator iset;
  	std::list<std::pair<int,int> > edgelist;
  	std::list<std::pair<int,int> >::iterator it;

	std::pair<int,int> edge;
	it = edgelist.begin();
	iset=edgeset.begin();

	if(mesh->elem2==NULL)
	    mesh->elem2=(int *)calloc(sizeof(int)*TETEDGE,mesh->ne);

	for (int eid=0;eid<mesh->ne;eid++)
	    for(int ed=0; ed<TETEDGE; ed++){
	        ee=(int*)(&mesh->elem[eid]);
		n1=MIN(ee[edgepair[ed][0]],ee[edgepair[ed][1]]);
		n2=MAX(ee[edgepair[ed][0]],ee[edgepair[ed][1]]);
		edge=std::make_pair(n1,n2);
		setlen=edgeset.size();
	        edgeset.insert(iset,edge);
		if(setlen<edgeset.size()){ // if no previous edge
		     edgelist.insert(it,edge);
		     pos=edgelist.size()-1;
		}else{
	             std::list<std::pair<int,int> >::iterator edidx = 
		         std::find(edgelist.begin(), edgelist.end(),edge);
		     pos = std::distance( edgelist.begin(), edidx) ;
		}
                mesh->elem2[eid*TETEDGE+ed]= pos;
	    }
	pos=0;

	mesh->nn+=edgelist.size();
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
