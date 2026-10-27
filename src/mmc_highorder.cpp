/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2026
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
#include <map>
#include <vector>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <string>
#include <string.h>
#include <sstream>

#include "mmc_mesh.h"
#include "mmc_highorder.h"
#include "mmc_utils.h"

#define TETEDGE 6

const int edgepair[TETEDGE][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

const int facelist[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

#ifdef __cplusplus
    extern "C"
#endif
void mesh_10nodetet(tetmesh* mesh, mcconfig* cfg) {
    int oldnn = mesh->nn;
    std::map< std::pair<int, int>, int > edgeidx;     /* edge (1-based node ids, low first) -> index of its mid-node */
    std::vector< std::pair<int, int> > edgelist;

    if (mesh->elem2 == NULL) {
        mesh->elem2 = (int*)calloc(sizeof(int) * TETEDGE, mesh->ne);
    }

    for (int eid = 0; eid < mesh->ne; eid++) {
        int* ee = mesh->elem + eid * mesh->elemlen;

        for (int ed = 0; ed < TETEDGE; ed++) {
            std::pair<int, int> edge = std::make_pair(MIN(ee[edgepair[ed][0]], ee[edgepair[ed][1]]), MAX(ee[edgepair[ed][0]], ee[edgepair[ed][1]]));
            std::pair<std::map< std::pair<int, int>, int >::iterator, bool> res = edgeidx.insert(std::make_pair(edge, (int)edgelist.size()));

            if (res.second) {
                edgelist.push_back(edge);
            }

            mesh->elem2[eid * TETEDGE + ed] = res.first->second;
        }
    }

    mesh->nn += edgelist.size();

    int issharednode = (cfg->node == mesh->node); /* a JSON mesh: cfg owns the node buffer, see mesh_clear */
    mesh->node = (FLOAT3*)realloc((void*)mesh->node, sizeof(FLOAT3) * (mesh->nn));

    if (issharednode) {
        cfg->node = mesh->node;
        cfg->nodenum = mesh->nn;
    }

    /* one slot per source or pattern, as mesh_loadelem allocates it */
    size_t wlen = (size_t)mesh->nn * cfg->maxgate * ((cfg->extrasrclen > cfg->srcnum) ? cfg->extrasrclen : cfg->srcnum);
    mesh->weight = (double*)realloc((void*)mesh->weight, sizeof(double) * wlen);
    memset(mesh->weight, 0, sizeof(double) * wlen);

    for (size_t k = 0; k < edgelist.size(); k++) {    /* the mid-node of each edge; node ids are 1-based */
        FLOAT3* a = mesh->node + edgelist[k].first - 1, *b = mesh->node + edgelist[k].second - 1, *m = mesh->node + oldnn + k;
        m->x = (a->x + b->x) * 0.5f;
        m->y = (a->y + b->y) * 0.5f;
        m->z = (a->z + b->z) * 0.5f;
    }
}


class Vec3Hash {
  public:
    std::size_t operator()(std::array<int, 3> const& vec) const {
        std::size_t seed = 3;

        for (int i = 0; i < 3; i++) {
            seed ^= vec[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        return seed;
    }
};

#ifdef __cplusplus
    extern "C"
#endif
void mesh_getfacenb(tetmesh* mesh, mcconfig* cfg) {
    std::unordered_map< std::array<int, 3>, std::pair<unsigned int, unsigned int>, Vec3Hash > facenb;
    std::unordered_map< std::array<int, 3>, std::pair<unsigned int, unsigned int>, Vec3Hash > ::iterator it;

    for (int i = 0; i < mesh->ne; i++) {
        int* ee = mesh->elem + i * mesh->elemlen;

        for (int j = 0; j < 4; j++) {
            std::array<int, 3> facevec = { ee[facelist[j][0]], ee[facelist[j][1]], ee[facelist[j][2]] };
            std::sort (facevec.begin(), facevec.end());

            if (facenb.find(facevec) != facenb.end()) {
                facenb[facevec].second = (i << 2) + j + 1;
            } else {
                facenb[facevec] = std::make_pair((i << 2) + j, 0);
            }
        }
    }

    if (mesh->facenb) {
        free(mesh->facenb);
    }

    mesh->facenb = (int*)calloc(sizeof(int) * mesh->elemlen, mesh->ne);

    for (it = facenb.begin(); it != facenb.end(); it++) {
        if (it->second.second != 0) {
            mesh->facenb[it->second.first] = ((it->second.second - 1) >> 2) + 1;
            mesh->facenb[it->second.second - 1] = (it->second.first >> 2) + 1;
        }
    }
}
