/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2021
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing
**          in Plucker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli,
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a>
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**  \li \c (\b Fang2019) Qianqian Fang and Shijie Yan,
**          <a href="http://dx.doi.org/10.1117/1.JBO.24.11.115002">
**          "Graphics processing unit-accelerated mesh-based Monte Carlo photon transport
**           simulations,"</a> J. of Biomedical Optics, 24(11), 115002 (2019)
**  \li \c (\b Yuan2021) Yaoshen Yuan, Shijie Yan, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-12-1-147">
**          "Light transport modeling in highly complex tissues using the implicit
**           mesh-based Monte Carlo algorithm,"</a> Biomed. Optics Express, 12(1) 147-161 (2021)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

#include <set>
#include <list>
#include <algorithm>
#include <string.h>
#include "mmc_mesh.h"
#include "mmc_highorder.h"
#include "mmc_utils.h"

#define TETEDGE 6

const int edgepair[TETEDGE][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

#ifdef __cplusplus
    extern "C"
#endif
void mesh_10nodetet(tetmesh* mesh, mcconfig* cfg) {
    int pos, n1, n2, oldnn = mesh->nn;
    unsigned int setlen;
    int* ee;

    std::set< std::pair<int, int> > edgeset;
    std::set< std::pair<int, int> >::iterator iset;
    std::list <std::pair<int, int> > edgelist;
    std::list <std::pair<int, int> >::iterator it;

    std::pair<int, int> edge;
    it = edgelist.begin();
    iset = edgeset.begin();

    if (mesh->elem2 == NULL) {
        mesh->elem2 = (int*)calloc(sizeof(int) * TETEDGE, mesh->ne);
    }

    for (int eid = 0; eid < mesh->ne; eid++)
        for (int ed = 0; ed < TETEDGE; ed++) {
            ee = (int*)(&mesh->elem[eid]);
            n1 = MIN(ee[edgepair[ed][0]], ee[edgepair[ed][1]]);
            n2 = MAX(ee[edgepair[ed][0]], ee[edgepair[ed][1]]);
            edge = std::make_pair(n1, n2);
            setlen = edgeset.size();
            edgeset.insert(iset, edge);

            if (setlen < edgeset.size()) { // if no previous edge
                edgelist.insert(it, edge);
                pos = edgelist.size() - 1;
            } else {
                std::list< std::pair<int, int> >::iterator edidx =
                    std::find(edgelist.begin(), edgelist.end(), edge);
                pos = std::distance( edgelist.begin(), edidx) ;
            }

            mesh->elem2[eid * TETEDGE + ed] = pos;
        }

    pos = 0;

    mesh->nn += edgelist.size();
    mesh->node = (MMCfloat3*)realloc((void*)mesh->node, sizeof(MMCfloat3) * (mesh->nn));
    mesh->weight = (double*)realloc((void*)mesh->weight, sizeof(double) * mesh->nn * cfg->maxgate);
    memset(mesh->weight, 0, sizeof(double)*mesh->nn * cfg->maxgate); // if mesh->weight is filled, need to allocate a new buffer, and copy the old buffer gate by gate

    for (it = edgelist.begin(); it != edgelist.end(); it++) {
        for (int i = 0; i < 3; i++) {
            ((float*)(&mesh->node[oldnn + pos]))[i] =
                (((float*)(&mesh->node[(*it).first]))[i] + ((float*)(&mesh->node[(*it).second]))[i]) * 0.5f;
        }

        pos++;
    }
}
