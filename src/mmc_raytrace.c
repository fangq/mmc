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

/***************************************************************************//**
\file    tettracing.c

\brief   Core unit for mash-based photon transport using ray tracing algorithms
*******************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mmc_raytrace.h"

/**<  Macro to enable SSE4 based ray-tracers */

#ifdef MMC_USE_SSE
    #include <smmintrin.h>
    __m128 int_coef;                            /**<  a global variable used for SSE4 ray-tracers */
#endif

#define F32N(a) ((a) & 0x80000000)          /**<  Macro to test if a floating point is negative */
#define F32P(a) ((a) ^ 0x80000000)          /**<  Macro to test if a floating point is positive */

/**
 * \mapping from edge index to node index
 */

const int e2n[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

/**
 * \brief Tetrahedron face edge indices
 *
 * fc[4] points to the 4 facets of a tetrahedron, with each
 * triangular face made of 3 directed edges. The numbers [0-5] are the
 * local indices of the edge, defined by
 * 0:[0->1], 1:[0->2], 2: [0->3], 3:[1->2], 4:[1->3], 5:[2->3]
 * where the pair in [] are the local nodes forming the edge
 *
 */

const int fc[4][3] = {{0, 4, 2}, {3, 5, 4}, {2, 5, 1}, {1, 3, 0}};

/**
 * \brief Tetrahedron faces, in clock-wise orders, represented using local node indices
 *
 * node-connectivity, i.e. nc[4] points to the 4 facets of a tetrahedron, with each
 * triangular face made of 3 nodes. The numbers [0-4] are the
 * local node indices (starting from 0). The order of the nodes
 * are in clock-wise orders.
 */

const int nc[4][3] = {{3, 0, 1}, {3, 1, 2}, {2, 0, 3}, {1, 0, 2}};

/**
 * \brief Tetrahedron faces, in counter-clock-wise orders, represented using local node indices
 *
 * out is like nc[] but with counter-clock-wise orientation. this makes
 * the normal vec of each face pointing outwards.
 */

extern const int out[4][3];     /**< defined in simpmesh.c, node orders for each face, in counter-clock-wise orders*/

/**
 * \brief The local index of the node with an opposite face to the i-th face defined in nc[][]
 *
 * nc[i] <-> node[facemap[i]]
 * the 1st face of this tet, i.e. nc[0]={3,0,1}, is opposite to the 3rd node
 * the 2nd face of this tet, i.e. nc[1]={3,1,2}, is opposite to the 1st node
 * etc.
 */

extern const int facemap[4];

/**
 * \brief Inverse mapping between the local index of the node and the corresponding opposite face in nc[]'s order
 *
 * nc[ifacemap[i]] <-> node[i]
 * the 1st node of this tet is in opposite to the 2nd face, i.e. nc[1]={3,1,2}
 * the 2nd node of this tet is in opposite to the 3rd face, i.e. nc[1]={2,0,3}
 * etc.
 */

extern const int ifacemap[4];

/**
 * \brief Index mapping from the i-th face-neighbors (facenb) to the face defined in nc[][]
 *
 * facenb[i] <-> nc[faceorder[i]]
 * the 1st tet neighbor shares the 2nd face of this tet, i.e. nc[1]={3,1,2}
 * the 2nd tet neighbor shares the 4th face of this tet, i.e. nc[3]={1,0,2}
 * etc.
 */

extern const int faceorder[5];

/**
 * \brief Index mapping from the i-th face defined in nc[][] to the face-neighbor (facenb) face orders
 *
 * nc[ifaceorder[i]] <-> facenb[i]
 * nc[0], made of nodes {3,0,1}, is the face connecting to the 4th neighbor (facenb[3]),
 * nc[1], made of nodes {3,1,2}, is the face connecting to the 1st neighbor (facenb[0]),
 * etc.
 */

extern const int ifaceorder[4];

/**
 * \brief A mapping from SSE4 mask output to the face index
 */

const char maskmap[16] = {4, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3};

/**
 * \brief function to linearly interpolate between 3 3D points (p1,p2,p3) using weight (w)
 *
 * pout=p1*w1+p2*w2+p3*w3
 * \param[in] w: weight vector
 * \param[in] p1: 1st point
 * \param[in] p2: 2nd point
 * \param[in] p3: 3rd point
 * \param[out] pout: output 3D position
 */

inline void interppos(MMCfloat3* w, MMCfloat3* p1, MMCfloat3* p2, MMCfloat3* p3, MMCfloat3* pout) {
    pout->x = w->x * p1->x + w->y * p2->x + w->z * p3->x;
    pout->y = w->x * p1->y + w->y * p2->y + w->z * p3->y;
    pout->z = w->x * p1->z + w->y * p2->z + w->z * p3->z;
}

/**
 * \brief function to linearly interpolate between 3 3D points (p1,p2,p3) using scalar weight (w)
 *
 * pout=p1*w1+p2*w2+p3*w3
 * \param[in] w1: weight for p1
 * \param[in] w2: weight for p2
 * \param[in] w3: weight for p3
 * \param[in] p1: 1st point
 * \param[in] p2: 2nd point
 * \param[in] p3: 3rd point
 * \param[out] pout: output 3D position
 */

inline void getinterp(float w1, float w2, float w3, MMCfloat3* p1, MMCfloat3* p2, MMCfloat3* p3, MMCfloat3* pout) {
    pout->x = w1 * p1->x + w2 * p2->x + w3 * p3->x;
    pout->y = w1 * p1->y + w2 * p2->y + w3 * p3->y;
    pout->z = w1 * p1->z + w2 * p2->z + w3 * p3->z;
}

/**
 * \brief Function to deal with ray-edge/ray-vertex intersections
 *
 * when a photon is crossing a vertex or edge, (slightly) pull the
 * photon toward the center of the element and try again
 *
 * \param[in,out] p: current photon position
 * \param[in] nodes: pointer to the 4 nodes of the tet
 * \param[in] ee: indices of the 4 nodes ee=elem[eid]
 */

void fixphoton(MMCfloat3* p, MMCfloat3* nodes, int* ee) {
    MMCfloat3 c0 = {0.f, 0.f, 0.f};
    int i;

    /*calculate element centroid*/
    for (i = 0; i < 4; i++) {
        vec_add(&c0, nodes + ee[i] - 1, &c0);
    }

    p->x += (c0.x * 0.25f - p->x) * FIX_PHOTON;
    p->y += (c0.y * 0.25f - p->y) * FIX_PHOTON;
    p->z += (c0.z * 0.25f - p->z) * FIX_PHOTON;
}

/**
 * \brief Plucker-coordinate based ray-tracer to advance photon by one step
 *
 * this function uses Plucker-coordinate based ray-triangle intersection
 * tests to advance photon by one step, see Fang2010.
 *
 * \param[in,out] r: the current ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

float plucker_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {

    MMCfloat3 pcrx = {0.f}, p1 = {0.f};
    medium* prop;
    int* ee;
    int i, tshift, eid, faceidx = -1;
    float w[6] = {0.f}, Rv, ww, currweight, dlen = 0.f, rc, mus; /*dlen is the physical distance*/
    unsigned int* wi = (unsigned int*)w;
    float baryout[4] = {0.f, 0.f, 0.f, 0.f}, *baryp0 = &(r->bary0.x);
    float Lp0 = 0.f, ratio;

    if (tracer->mesh == NULL || tracer->d == NULL || r->eid <= 0 || r->eid > tracer->mesh->ne) {
        return -1;
    }

    eid = r->eid - 1;
    r->faceid = -1;
    r->isend = 0;
    r->Lmove = 0.f;
    vec_add(&(r->p0), &(r->vec), &p1);
    vec_cross(&(r->p0), &p1, &pcrx);
    ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);
    prop = tracer->mesh->med + (tracer->mesh->type[eid]);
    rc = prop->n * R_C0;
    currweight = r->weight;
    mus = (cfg->mcmethod == mmMCX) ? prop->mus : (prop->mua + prop->mus);

#ifdef MMC_USE_SSE
    {
        __m128 O = _mm_load_ps(&(r->vec.x));
        __m128 M, D, T = _mm_load_ps(&(pcrx.x));

        for (i = 0; i < 6; i++) {
            D = _mm_load_ps(&(tracer->d[(eid) * 6 + i].x));
            M = _mm_load_ps(&(tracer->m[(eid) * 6 + i].x));
#ifdef __SSE4_1__
            w[i] = _mm_cvtss_f32(_mm_add_ss(_mm_dp_ps(O, M, 0x7F), _mm_dp_ps(T, D, 0x7F)));
#else
            M = _mm_mul_ps(O, M);
            M = _mm_hadd_ps(M, M);
            M = _mm_hadd_ps(M, M);
            D = _mm_mul_ps(T, D);
            D = _mm_hadd_ps(D, D);
            D = _mm_hadd_ps(D, D);
            _mm_store_ss(w + i, _mm_add_ss(M, D));
#endif
        }
    }
#else

    for (i = 0; i < 6; i++) {
        w[i] = pinner(&(r->vec), &pcrx, tracer->d + (eid) * 6 + i, tracer->m + (eid) * 6 + i);
    }

#endif

    //exit(1);
    if (cfg->debuglevel & dlTracing) {
        MMC_FPRINTF(cfg->flog, "%d \n", eid);
    }

    r->pout.x = MMC_UNDEFINED;

    for (i = 0; i < 4; i++) {
        //if(cfg->debuglevel&dlTracing) MMC_FPRINTF(cfg->flog,"testing face [%d]\n",i);

        if (i >= 2) {
            w[fc[i][1]] = -w[fc[i][1]];    // can not go back
        }

        if (F32N(wi[fc[i][0]]) & F32N(wi[fc[i][1]]) & F32P(wi[fc[i][2]])) {
            // f_leave
            //if(cfg->debuglevel&dlTracingExit) MMC_FPRINTF(cfg->flog,"ray exits face %d[%d] of %d\n",i,faceorder[i],eid);

            Rv = 1.f / (-w[fc[i][0]] - w[fc[i][1]] + w[fc[i][2]]);
            baryout[nc[i][0]] = -w[fc[i][0]] * Rv;
            baryout[nc[i][1]] = -w[fc[i][1]] * Rv;
            baryout[nc[i][2]] = w[fc[i][2]] * Rv;

            getinterp(baryout[nc[i][0]], baryout[nc[i][1]], baryout[nc[i][2]],
                      tracer->mesh->node + ee[nc[i][0]] - 1,
                      tracer->mesh->node + ee[nc[i][1]] - 1,
                      tracer->mesh->node + ee[nc[i][2]] - 1, &(r->pout));

            //if(cfg->debuglevel&dlTracingExit) MMC_FPRINTF(cfg->flog,"exit point %f %f %f\n",r->pout.x,r->pout.y,r->pout.z);

            Lp0 = dist(&(r->p0), &(r->pout));
            dlen = (mus <= EPS) ? R_MIN_MUS : r->slen / mus;
            faceidx = i;
            r->faceid = faceorder[i];

            {
                int* enb = (int*)(tracer->mesh->facenb + eid * tracer->mesh->elemlen);
                r->nexteid = enb[r->faceid];
#ifdef MMC_USE_SSE
                int nexteid = (r->nexteid - 1) * 6;

                if (r->nexteid > 0) {
                    _mm_prefetch((char*) & ((tracer->m + nexteid)->x), _MM_HINT_T0);
                    _mm_prefetch((char*) & ((tracer->m + (nexteid) * 6 + 4)->x), _MM_HINT_T0);
                    _mm_prefetch((char*) & ((tracer->d + (nexteid) * 6)->x), _MM_HINT_T0);
                    _mm_prefetch((char*) & ((tracer->d + (nexteid) * 6 + 4)->x), _MM_HINT_T0);
                }

#endif
            }
            r->isend = (Lp0 > dlen);
            r->Lmove = ((r->isend) ? dlen : Lp0);
            break;
        }
    }

    visit->raytet++;

    if (tracer->mesh->type[eid] == 0) {
        visit->raytet0++;
    }

    if (r->pout.x != MMC_UNDEFINED) {
        if ((int)((r->photontimer + r->Lmove * rc - cfg->tstart)*visit->rtstep) >= (int)((cfg->tend - cfg->tstart)*visit->rtstep)) { /*exit time window*/
            r->faceid = -2;
            r->pout.x = MMC_UNDEFINED;
            r->Lmove = (cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f;
        }

        if (cfg->mcmethod == mmMCX) {
            r->weight *= expf(-prop->mua * r->Lmove);
        }

        if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otJacobian) {
            currweight = expf(-DELTA_MUA * r->Lmove);
            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        } else if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWL) {
            currweight = r->Lmove;
            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        }

        r->slen -= r->Lmove * mus;

        if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWP) {
            if (r->slen0 < EPS) {
                currweight = 1;
            } else {
                currweight = r->Lmove * mus / r->slen0;
            }

            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        }

        r->p0.x += r->Lmove * r->vec.x;
        r->p0.y += r->Lmove * r->vec.y;
        r->p0.z += r->Lmove * r->vec.z;
        //if(cfg->debuglevel&dlWeight) MMC_FPRINTF(cfg->flog,"update weight to %f and path end %d \n",r->weight,r->isend);

        if (!cfg->basisorder) {
            ww = currweight - r->weight;
            r->Eabsorb += ww;
            r->photontimer += r->Lmove * rc;

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * tracer->mesh->ne;
            } else {
                tshift = MIN( ((int)((r->photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * tracer->mesh->ne;
            }

            if (cfg->mcmethod == mmMCX) {
                if (cfg->srctype != stPattern) {
                    if (cfg->isatomic)
                        #pragma omp atomic
                        tracer->mesh->weight[eid + tshift] += ww;
                    else {
                        tracer->mesh->weight[eid + tshift] += ww;
                    }
                } else {
                    int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                    int pidx; // pattern index

                    for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                        if (cfg->isatomic)
                            #pragma omp atomic
                            tracer->mesh->weight[(eid + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx];
                        else {
                            tracer->mesh->weight[(eid + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx];
                        }
                    }
                }
            }
        } else {
            if (cfg->debuglevel & dlBary) MMC_FPRINTF(cfg->flog, "Y [%f %f %f %f]\n",
                        baryout[0], baryout[1], baryout[2], baryout[3]);

            if (Lp0 > EPS) {
                r->photontimer += r->Lmove * rc;
                ww = currweight - r->weight;

                if (prop->mua > 0.f) {
                    ratio = r->Lmove / Lp0;
                    r->Eabsorb += ww;

                    if (cfg->outputtype != otEnergy && cfg->outputtype != otWP) {
                        ww /= prop->mua;
                    }

                    if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                        tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * tracer->mesh->nn;
                    } else {
                        tshift = MIN( ((int)((r->photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * tracer->mesh->nn;
                    }

                    if (cfg->debuglevel & dlAccum) MMC_FPRINTF(cfg->flog, "A %f %f %f %e %d %e\n",
                                r->p0.x - (r->Lmove * 0.5f)*r->vec.x, r->p0.y - (r->Lmove * 0.5f)*r->vec.y, r->p0.z - (r->Lmove * 0.5f)*r->vec.z, ww, eid + 1, dlen);

                    ww *= 0.5f;

                    if (r->isend)
                        for (i = 0; i < 4; i++) {
                            baryout[i] = (1.f - ratio) * baryp0[i] + ratio * baryout[i];
                        }

                    if (cfg->mcmethod == mmMCX) {
                        if (cfg->srctype != stPattern) {
                            if (cfg->isatomic)
                                for (i = 0; i < 4; i++)
                                    #pragma omp atomic
                                    tracer->mesh->weight[ee[i] - 1 + tshift] += ww * (baryp0[i] + baryout[i]);
                            else
                                for (i = 0; i < 4; i++) {
                                    tracer->mesh->weight[ee[i] - 1 + tshift] += ww * (baryp0[i] + baryout[i]);
                                }
                        } else {
                            int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                            int pidx; // pattern index

                            for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                if (cfg->isatomic)
                                    for (i = 0; i < 4; i++)
                                        #pragma omp atomic
                                        tracer->mesh->weight[(ee[i] - 1 + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx] * (baryp0[i] + baryout[i]);
                                else
                                    for (i = 0; i < 4; i++) {
                                        tracer->mesh->weight[(ee[i] - 1 + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx] * (baryp0[i] + baryout[i]);
                                    }
                            }
                        }
                    }
                }

                if (r->isend) {
                    memcpy(baryp0, baryout, sizeof(float4));
                } else {
                    if (r->nexteid && faceidx >= 0) {
                        int j, k, *nextenb = (int*)(tracer->mesh->elem + (r->nexteid - 1) * tracer->mesh->elemlen);
                        memset(baryp0, 0, sizeof(float4));

                        for (j = 0; j < 3; j++)
                            for (k = 0; k < 4; k++) {
                                if (ee[nc[faceidx][j]] == nextenb[k]) {
                                    baryp0[k] = baryout[nc[faceidx][j]];
                                    break;
                                }
                            }
                    }
                }
            }
        }

        if (r->faceid == -2) {
            return 0.f;
        }
    }

    return r->slen;
}

#ifdef MMC_USE_SSE

inline __m128 rcp_nr(const __m128 a) {
    const __m128 r = _mm_rcp_ps(a);
    return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, a), r));
}

#if !defined(__EMSCRIPTEN__)

/**
 * \brief Havel-based SSE4 ray-triangle intersection test
 *
 * this function uses Havel-based algorithm to test if a ray intersects
 * with a triangle see Fang2012.
 *
 * \param[in] vecN: the normal of the face
 * \param[out] bary: output the Barycentric coordinates of the intersection point
 * \param[in] o: current ray origin
 * \param[out] d: current ray direction
 */

inline int havel_sse4(MMCfloat3* vecN, MMCfloat3* bary, const __m128 o, const __m128 d) {
    const __m128 n = _mm_load_ps(&vecN->x);
    const __m128 det = _mm_dp_ps(n, d, 0x7f);

    if (_mm_movemask_ps(det) & 1) {
        return 0;
    }

    const __m128 dett = _mm_dp_ps(_mm_mul_ps(int_coef, n), o, 0xff);
    const __m128 oldt = _mm_load_ss(&bary->x);

    if ((_mm_movemask_ps(_mm_xor_ps(dett, _mm_sub_ss(_mm_mul_ss(oldt, det), dett))) & 1) == 0) {
        const __m128 detp = _mm_add_ps(_mm_mul_ps(o, det), _mm_mul_ps(dett, d));
        const __m128 detu = _mm_dp_ps(detp, _mm_load_ps(&((vecN + 1)->x)), 0xf1);

        if ((_mm_movemask_ps(_mm_xor_ps(detu, _mm_sub_ss(det, detu))) & 1) == 0) {
            const __m128 detv = _mm_dp_ps(detp, _mm_load_ps(&((vecN + 2)->x)), 0xf1);

            if ((_mm_movemask_ps(_mm_xor_ps(detv, _mm_sub_ss(det, _mm_add_ss(detu, detv)))) & 1) == 0) {
                const __m128 inv_det = rcp_nr(det);
                _mm_store_ss(&bary->x, _mm_mul_ss(dett, inv_det));
                _mm_store_ss(&bary->y, _mm_mul_ss(detu, inv_det));
                _mm_store_ss(&bary->z, _mm_mul_ss(detv, inv_det));
                //_mm_store_ps(&pout->x, _mm_mul_ps(detp,_mm_shuffle_ps(inv_det, inv_det, 0)));
                return (bary->x == bary->x); // when a photon is outside, bary->x=NaN
            }
        }
    }

    return 0;
}

/**
 * \brief Havel-based SSE4 ray-tracer to advance photon by one step
 *
 * this function uses Havel-based SSE4 ray-triangle intersection
 * tests to advance photon by one step, see Fang2012.
 *
 * \param[in,out] r: the current ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

float havel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {

    MMCfloat3 bary = {1e10f, 0.f, 0.f, 0.f};
    float barypout[4] __attribute__ ((aligned(16)));
    float rc, currweight, dlen, ww, Lp0, mus;
    int i, j, k, tshift, *enb = NULL, *nextenb = NULL, eid;
    __m128 O, T, S;

    if (tracer->mesh == NULL || tracer->m == NULL || r->eid <= 0 || r->eid > tracer->mesh->ne) {
        return -1;
    }

    r->p0.w = 1.f;
    r->vec.w = 0.f;
    eid = r->eid - 1;

    O = _mm_load_ps(&(r->p0.x));
    T = _mm_load_ps(&(r->vec.x));

    medium* prop;
    int* ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);
    prop = tracer->mesh->med + (tracer->mesh->type[eid]);
    rc = prop->n * R_C0;
    currweight = r->weight;
    mus = (cfg->mcmethod == mmMCX) ? prop->mus : (prop->mua + prop->mus);

    r->pout.x = MMC_UNDEFINED;
    r->faceid = -1;
    r->isend = 0;
    r->Lmove = 0.f;

    for (i = 0; i < 4; i++)
        if (havel_sse4(tracer->m + eid * 12 + i * 3, &bary, O, T)) {

            r->faceid = faceorder[i];
            dlen = (mus <= EPS) ? R_MIN_MUS : r->slen / mus;
            Lp0 = bary.x;
            r->isend = (Lp0 > dlen);
            r->Lmove = ((r->isend) ? dlen : Lp0);

            if (!r->isend) {
                enb = (int*)(tracer->mesh->facenb + eid * tracer->mesh->elemlen);
                r->nexteid = enb[r->faceid];

                if (r->nexteid > 0) {
                    nextenb = (int*)(tracer->m + (r->nexteid - 1) * 12);
                    _mm_prefetch((char*)(nextenb), _MM_HINT_T0);
                    _mm_prefetch((char*)(nextenb + 4 * 16), _MM_HINT_T0);
                    _mm_prefetch((char*)(nextenb + 8 * 16), _MM_HINT_T0);
                    nextenb = (int*)(tracer->mesh->elem + (r->nexteid - 1) * tracer->mesh->elemlen);
                }
            }

            S = _mm_set1_ps(bary.x);
            S = _mm_mul_ps(S, T);
            S = _mm_add_ps(S, O);
            _mm_store_ps(&(r->pout.x), S);

            if ((int)((r->photontimer + r->Lmove * rc - cfg->tstart)*visit->rtstep) >= (int)((cfg->tend - cfg->tstart)*visit->rtstep)) { /*exit time window*/
                r->faceid = -2;
                r->pout.x = MMC_UNDEFINED;
                r->Lmove = (cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f;
            }

            if (cfg->mcmethod == mmMCX) {
                r->weight *= expf(-prop->mua * r->Lmove);
            }

            if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otJacobian) {
                currweight = expf(-DELTA_MUA * r->Lmove);
                currweight *= cfg->replayweight[r->photonid];
                currweight += r->weight;
            } else if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWL) {
                currweight = r->Lmove;
                currweight *= cfg->replayweight[r->photonid];
                currweight += r->weight;
            }

            r->slen -= r->Lmove * mus;

            if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWP) {
                if (r->slen0 < EPS) {
                    currweight = 1;
                } else {
                    currweight = r->Lmove * mus / r->slen0;
                }

                currweight *= cfg->replayweight[r->photonid];
                currweight += r->weight;
            }

            if (bary.x == 0.f) {
                break;
            }

            ww = currweight - r->weight;
            r->photontimer += r->Lmove * rc;

            if (prop->mua > 0.f) {
                r->Eabsorb += ww;

                if (cfg->outputtype != otEnergy && cfg->outputtype != otWP) {
                    ww /= prop->mua;
                }
            }

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            } else {
                tshift = MIN( ((int)((r->photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            }

            if (cfg->debuglevel & dlAccum) MMC_FPRINTF(cfg->flog, "A %f %f %f %e %d %e\n",
                        r->p0.x, r->p0.y, r->p0.z, bary.x, eid + 1, dlen);

            T = _mm_mul_ps(T, _mm_set1_ps(r->Lmove));
            S = _mm_add_ps(O, T);
            _mm_store_ps(&(r->p0.x), S);

            barypout[out[i][0]] = 1.f - bary.y - bary.z;
            barypout[out[i][1]] = bary.y;
            barypout[out[i][2]] = bary.z;
            barypout[facemap[i]] = 0.f;

            if (cfg->debuglevel & dlBary) MMC_FPRINTF(cfg->flog, "Y [%f %f %f %f]\n",
                        barypout[0], barypout[1], barypout[2], barypout[3]);

            T = _mm_load_ps(barypout);      /* bary centric at pout */
            O = _mm_load_ps(&(r->bary0.x)); /* bary centric at p0 */

            dlen = r->Lmove / bary.x;       /* normalized moving length */

            if (cfg->debuglevel & dlAccum) MMC_FPRINTF(cfg->flog, "A %f %f %f %e %d %e\n",
                        r->p0.x - (r->Lmove * 0.5f)*r->vec.x, r->p0.y - (r->Lmove * 0.5f)*r->vec.y, r->p0.z - (r->Lmove * 0.5f)*r->vec.z, ww, eid + 1, dlen);

            if (r->isend) {                 /* S is the bary centric for the photon after move */
                S = _mm_add_ps(_mm_mul_ps(T, _mm_set1_ps(dlen)), _mm_mul_ps(O, _mm_set1_ps(1.f - dlen)));
            } else {
                S = T;
            }

            //if(cfg->debuglevel&dlBary)
            //    MMC_FPRINTF(cfg->flog,"old bary0=[%f %f %f %f]\n",r->bary0.x,r->bary0.y,r->bary0.z,r->bary0.w);

            _mm_store_ps(barypout, S);

            if (nextenb && enb) {
                float* barynext = &(r->bary0.x);
                _mm_store_ps(barynext, _mm_set1_ps(0.f));

                for (j = 0; j < 3; j++)
                    for (k = 0; k < 4; k++) {
                        if (ee[out[i][j]] == nextenb[k]) {
                            barynext[k] = barypout[out[i][j]];
                            break;
                        }
                    }

                //if(cfg->debuglevel&dlBary) MMC_FPRINTF(cfg->flog,"[%d %d %d %d],[%d %d %d %d] - ",
                //       ee[0],ee[1],ee[2],ee[3],nextenb[0],nextenb[1],nextenb[2],nextenb[3]);
                //if(cfg->debuglevel&dlBary) MMC_FPRINTF(cfg->flog,"[%f %f %f %f],[%f %f %f %f]\n",barypout[0],barypout[1],barypout[2],barypout[3],
                //       barynext[0],barynext[1],barynext[2],barynext[3]);
            } else {
                _mm_store_ps(&(r->bary0.x), S);
            }

            //if(cfg->debuglevel&dlBary)
            //    MMC_FPRINTF(cfg->flog,"new bary0=[%f %f %f %f]\n",r->bary0.x,r->bary0.y,r->bary0.z,r->bary0.w);
            if (cfg->mcmethod == mmMCX) {
                if (!cfg->basisorder) {
                    if (cfg->srctype != stPattern) {
                        if (cfg->isatomic)
                            #pragma omp atomic
                            tracer->mesh->weight[eid + tshift] += ww;
                        else {
                            tracer->mesh->weight[eid + tshift] += ww;
                        }
                    } else {
                        int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                        int pidx; // pattern index

                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                            if (cfg->isatomic)
                                #pragma omp atomic
                                tracer->mesh->weight[(eid + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx];
                            else {
                                tracer->mesh->weight[(eid + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx];
                            }
                        }
                    }
                } else {
                    T = _mm_mul_ps(_mm_add_ps(O, S), _mm_set1_ps(ww * 0.5f));
                    _mm_store_ps(barypout, T);

                    if (cfg->srctype != stPattern) {
                        if (cfg->isatomic)
                            for (j = 0; j < 4; j++)
                                #pragma omp atomic
                                tracer->mesh->weight[ee[j] - 1 + tshift] += barypout[j];
                        else
                            for (j = 0; j < 4; j++) {
                                tracer->mesh->weight[ee[j] - 1 + tshift] += barypout[j];
                            }
                    } else {
                        int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                        int pidx; // pattern index

                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                            if (cfg->isatomic)
                                for (j = 0; j < 4; j++)
                                    #pragma omp atomic
                                    tracer->mesh->weight[(ee[j] - 1 + tshift)*cfg->srcnum + pidx] += barypout[j] * cfg->srcpattern[pidx * psize + r->posidx];
                            else
                                for (j = 0; j < 4; j++) {
                                    tracer->mesh->weight[(ee[j] - 1 + tshift)*cfg->srcnum + pidx] += barypout[j] * cfg->srcpattern[pidx * psize + r->posidx];
                                }
                        }
                    }
                }
            }

            break;
        }

    visit->raytet++;

    if (tracer->mesh->type[eid] == 0) {
        visit->raytet0++;
    }

    r->p0.w = 0.f;

    if (r->faceid == -2) {
        return 0.f;
    }

    return r->slen;
}

/**
 * \brief Badouel-based SSE4 ray-tracer to advance photon by one step
 *
 * this function uses Badouel-based SSE4 ray-triangle intersection
 * tests to advance photon by one step, see Fang2012. Both Badouel and
 * Branch-less Badouel algorithms do not calculate the Barycentric coordinates
 * and can only store energy loss using 0-th order basis function.
 *
 * \param[in,out] r: the current ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

float badouel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {

    MMCfloat3 bary = {1e10f, 0.f, 0.f, 0.f};
    float Lp0 = 0.f, rc, currweight, dlen, ww, t[4] = {1e10f, 1e10f, 1e10f, 1e10f};
    int i, tshift, faceidx = -1, eid;

    if (tracer->mesh == NULL || tracer->m == NULL || r->eid <= 0 || r->eid > tracer->mesh->ne) {
        return -1;
    }

    r->p0.w = 1.f;
    r->vec.w = 0.f;
    eid = r->eid - 1;

    r->pout.x = MMC_UNDEFINED;
    r->faceid = -1;
    r->isend = 0;

    const __m128 o = _mm_load_ps(&(r->p0.x));
    const __m128 d = _mm_load_ps(&(r->vec.x));

    for (i = 0; i < 4; i++) {
        __m128 det, dett, n;
        __m128 inv_det;
        n = _mm_load_ps(&(tracer->m[3 * ((eid << 2) + i)].x));
        det = _mm_dp_ps(n, d, 0x7f);

        if (_mm_movemask_ps(det) & 1) {
            continue;
        }

        dett = _mm_dp_ps(_mm_mul_ps(int_coef, n), o, 0xff);
        inv_det = rcp_nr(det);
        _mm_store_ss(t + i, _mm_mul_ss(dett, inv_det));

        if (t[i] < bary.x) {
            bary.x = t[i];
            faceidx = i;
            r->faceid = faceorder[i];
        }
    }

    if (r->faceid >= 0) {
        medium* prop;
        int* ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);
        float mus;
        prop = tracer->mesh->med + (tracer->mesh->type[eid]);
        rc = prop->n * R_C0;
        currweight = r->weight;

        mus = (cfg->mcmethod == mmMCX) ? prop->mus : (prop->mua + prop->mus);

        dlen = (mus <= EPS) ? R_MIN_MUS : r->slen / mus;
        Lp0 = bary.x;
        r->isend = (Lp0 > dlen);
        r->Lmove = ((r->isend) ? dlen : Lp0);

        r->pout.x = r->p0.x + bary.x * r->vec.x;
        r->pout.y = r->p0.y + bary.x * r->vec.y;
        r->pout.z = r->p0.z + bary.x * r->vec.z;

        if ((int)((r->photontimer + r->Lmove * rc - cfg->tstart)*visit->rtstep) >= (int)((cfg->tend - cfg->tstart)*visit->rtstep)) { /*exit time window*/
            r->faceid = -2;
            r->pout.x = MMC_UNDEFINED;
            r->Lmove = (cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f;
        }

        if (cfg->mcmethod == mmMCX) {
            r->weight *= expf(-prop->mua * r->Lmove);
        }

        if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otJacobian) {
            currweight = expf(-DELTA_MUA * r->Lmove);
            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        } else if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWL) {
            currweight = r->Lmove;
            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        }

        r->slen -= r->Lmove * mus;

        if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWP) {
            if (r->slen0 < EPS) {
                currweight = 1;
            } else {
                currweight = r->Lmove * mus / r->slen0;
            }

            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        }

        if (bary.x >= 0.f) {
            ww = currweight - r->weight;
            r->Eabsorb += ww;
            r->photontimer += r->Lmove * rc;

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            } else {
                tshift = MIN( ((int)((r->photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            }

            if (cfg->debuglevel & dlAccum) MMC_FPRINTF(cfg->flog, "A %f %f %f %e %d %e\n",
                        r->p0.x, r->p0.y, r->p0.z, bary.x, eid + 1, dlen);

            r->p0.x += r->Lmove * r->vec.x;
            r->p0.y += r->Lmove * r->vec.y;
            r->p0.z += r->Lmove * r->vec.z;

            if (cfg->mcmethod == mmMCX) {
                if (!cfg->basisorder) {
                    if (cfg->isatomic)
                        #pragma omp atomic
                        tracer->mesh->weight[eid + tshift] += ww;
                    else {
                        tracer->mesh->weight[eid + tshift] += ww;
                    }
                } else {
                    if (cfg->outputtype != otEnergy && cfg->outputtype != otWP) {
                        ww /= prop->mua;
                    }

                    ww *= 1.f / 3.f;

                    if (cfg->isatomic)
                        for (i = 0; i < 3; i++)
                            #pragma omp atomic
                            tracer->mesh->weight[ee[out[faceidx][i]] - 1 + tshift] += ww;
                    else
                        for (i = 0; i < 3; i++) {
                            tracer->mesh->weight[ee[out[faceidx][i]] - 1 + tshift] += ww;
                        }
                }
            }
        }
    }

    visit->raytet++;

    if (tracer->mesh->type[eid] == 0) {
        visit->raytet0++;
    }

    r->p0.w = 0.f;

    if (r->faceid == -2) {
        return 0.f;
    }

    return r->slen;
}

#else
float badouel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {
    MMC_ERROR(-6, "wrong option, please recompile with SSE4 enabled");
    return MMC_UNDEFINED;
}
float havel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {
    MMC_ERROR(-6, "wrong option, please recompile with SSE4 enabled");
    return MMC_UNDEFINED;
}
#endif  /* #if !defined(__EMSCRIPTEN__) */


/**
 * \brief Compute distances to edge (cylindrical) ROIs in the element
 *
 * Commpute the distance between P0 and the labeled edge
 * and the distance between P1 and the labeled edge.
 * P0 is the entry point and P1 is the exit point
 */
void compute_distances_to_edge(ray* r, raytracer* tracer, int* ee, int edgeid, float d2d[2], MMCfloat3 p2d[2], int* hitstatus) {
    MMCfloat3 u, OP;
    float r2;

    p2d[1] = tracer->mesh->node[ee[e2n[edgeid][0]] - 1];
    p2d[0] = tracer->mesh->node[ee[e2n[edgeid][1]] - 1];

    vec_diff(p2d + 1, p2d, &u);
    r2 = dist(p2d + 1, p2d); // get coordinates and compute distance between two nodes
    r2 = 1.0f / r2;
    vec_mult(&u, r2, &u); // normalized vector of the edge

    vec_diff(p2d + 1, &r->p0, &OP);
    r2 = vec_dot(&OP, &u);
    vec_mult(&u, r2, p2d);
    vec_diff(p2d, &OP, p2d);    // PP0: P0 projection on plane
    d2d[0] = vec_dot(p2d, p2d);

    vec_mult(&r->vec, r->Lmove, &OP);
    vec_add(&r->p0, &OP, &OP);      // P1
    vec_diff(p2d + 1, &OP, &OP);
    r2 = vec_dot(&OP, &u);
    vec_mult(&u, r2, p2d + 1);
    vec_diff(p2d + 1, &OP, p2d + 1); // PP1: P1 projection on plane
    d2d[1] = vec_dot(p2d + 1, p2d + 1);

    r2 = r->roisize[edgeid] * r->roisize[edgeid];
    *hitstatus = htNone;

    if (d2d[0] > r2 + EPS2 && d2d[1] < r2 - EPS2) {
        *hitstatus = htOutIn;
    } else if (d2d[0] < r2 - EPS2 && d2d[1] > r2 + EPS2) {
        *hitstatus = htInOut;
    } else if (d2d[1] > r2 + EPS2) {
        *hitstatus = htNoHitOut;
    } else if (d2d[1] < r2 - EPS2) {
        *hitstatus = htNoHitIn;
    }
}

/**
 * \brief Compute the next step length for edge-based iMMC
 */
float ray_cylinder_intersect(ray* r, int index, float d2d[2], MMCfloat3 p2d[2], int hitstatus) {
    float pt0[2], pt1[2], ph0[2], ph1[2], Lratio, theta;

    // update Lmove and reflection
    float dx, dy, dr2, delta, rt, sgn, Dp;
    rt = r->roisize[index];

    // accurate intersection computation
    if (hitstatus == htOutIn || hitstatus == htInOut) { // hit edgeroi out->in or in->out
        d2d[0] = sqrtf(d2d[0]);
        d2d[1] = sqrtf(d2d[1]);

        theta = vec_dot(p2d, p2d + 1);

        theta = theta / (d2d[0] * d2d[1] + EPS);
        pt0[0] = d2d[0];
        pt0[1] = 0.f;
        pt1[0] = d2d[1] * theta;
        pt1[1] = d2d[1] * sqrtf(fabs(1 - theta * theta));

        dx = pt1[0] - pt0[0];
        dy = pt1[1] - pt0[1];
        dr2 = dx * dx + dy * dy;
        Dp = pt0[0] * pt1[1];
        delta = sqrtf(rt * rt * dr2 - Dp * Dp); // must be >0
        dr2 = 1.0f / dr2;
        sgn = (dy >= 0) ? 1.0f : -1.0f;
        ph0[0] = (Dp * dy  + sgn * dx * delta) * dr2;
        ph0[1] = (-Dp * dx + fabs(dy) * delta) * dr2;
        ph1[0] = (Dp * dy  - sgn * dx * delta) * dr2;
        ph1[1] = (-Dp * dx - fabs(dy) * delta) * dr2;
        Dp = dist2d2(pt0, pt1);

        if (d2d[0] > rt && d2d[1] < rt) { // out->in
            Lratio = dist2d2(pt0, ph1);

            if (Dp > Lratio) {
                Lratio = sqrtf(Lratio / Dp);
            } else {
                Lratio = dist2d2(pt0, ph0);
                Lratio = sqrtf(Lratio / Dp);
            }
        } else { // d2d[0]<rt && d2d[1]>rt, in->out
            Lratio = dist2d2(pt1, ph1);

            if (Dp > Lratio) {
                Lratio = 1.f - sqrtf(Lratio / Dp);
            } else {
                Lratio = dist2d2(pt1, ph0);
                Lratio = 1.f - sqrtf(Lratio / Dp);
            }
        }

        return Lratio;
    }

    return 1.f;
}

/**
 * \brief Compute distances to node (spherical) ROIs if defined in the element
 *
 * Compute the distance between P0 and the labeled node
 * and the distance between P1 and the labeled node.
 * P0 is the entry point and P1 is the exit point
 */

void compute_distances_to_node(ray* r, raytracer* tracer, int* ee, int index, float nr, MMCfloat3** center, int* hitstatus) {
    float npdist0, npdist1;
    MMCfloat3 PP;

    nr = nr * nr;
    *center = tracer->mesh->node + ee[index] - 1;

    npdist0 = dist2(&r->p0, *center);   // P0
    vec_mult(&r->vec, r->Lmove, &PP);
    vec_add(&r->p0, &PP, &PP);
    npdist1 = dist2(&PP, *center);      // P1

    if (npdist0 > nr + EPS2 && npdist1 < nr - EPS2) {
        *hitstatus = htOutIn;
    } else if (npdist0 < nr - EPS2 && npdist1 > nr + EPS2) {
        *hitstatus = htInOut;
    } else if (npdist1 > nr + EPS2) {
        *hitstatus = htNoHitOut;
    } else if (npdist1 < nr - EPS2) {
        *hitstatus = htNoHitIn;
    } else {
        *hitstatus = htNone;
    }
}

/**
 * \brief Compute the next step length for node-based iMMC
 */
float ray_sphere_intersect(ray* r, int index, MMCfloat3* center, float nr, int hitstatus) {
    float temp1, temp2, d1, d2;
    MMCfloat3 oc;

    if (hitstatus == htOutIn || hitstatus == htInOut) {
        vec_diff(center, &r->p0, &oc);
        temp1 = vec_dot(&r->vec, &oc);
        temp2 = oc.x * oc.x + oc.y * oc.y + oc.z * oc.z - (nr * nr);
        temp2 = sqrtf(fabs(temp1 * temp1 - temp2));
        d1 = -temp1 + temp2;
        d2 = -temp1 - temp2;
        return ((hitstatus == htOutIn) ? MIN(d1, d2) : MAX(d1, d2));
    }

    return r->Lmove;
}

/**
 * \brief Compute the next step length for the face-based iMMC
 */
float ray_face_intersect(ray* r, raytracer* tracer, int* ee, int faceid, int baseid, int eid, int* hitstatus) {
    MMCfloat3 pf1, pv, fnorm, ptemp;
    float distf0, distf1, thick, ratio = 1.f;

    thick = r->roisize[faceid];

    vec_mult(&r->vec, r->Lmove, &ptemp);
    vec_add(&r->p0, &ptemp, &ptemp);    // P1: ptemp

    pf1 = tracer->mesh->node[ee[nc[ifaceorder[faceid]][0]] - 1]; // any point on face
    fnorm.x = (&(tracer->n[baseid].x))[ifaceorder[faceid]];     // normal vector of the face
    fnorm.y = (&(tracer->n[baseid].x))[ifaceorder[faceid] + 4];
    fnorm.z = (&(tracer->n[baseid].x))[ifaceorder[faceid] + 8];

    vec_diff(&r->p0, &pf1, &pv);
    distf0 = vec_dot(&pv, &fnorm);

    vec_diff(&ptemp, &pf1, &pv);
    distf1 = vec_dot(&pv, &fnorm);

    *hitstatus = htNone;

    if (distf0 > thick + EPS2 && distf1 < thick - EPS2) { // hit: out->in
        ratio = 1.f - (thick - distf1) / (distf0 - distf1);
        *hitstatus = htOutIn;
    } else if (distf0 < thick - EPS2 && distf1 > thick + EPS2) { // hit: in->out
        ratio = (thick - distf0) / (distf1 - distf0);
        *hitstatus = htInOut;
    } else if ((distf0 <= thick && distf1 >= thick) || (distf0 > thick && distf1 > thick)) {
        *hitstatus = htNoHitOut;
    } else if ((distf0 >= thick && distf1 <= thick) || (distf0 < thick && distf1 < thick)) {
        *hitstatus = htNoHitIn;
    }

    return ratio;
}

/**
 * \brief Implicit MMC ray-tracing core function
 *
 * if doinit is set to 1, this function only initialize r->inroi for 3 roi types
 * if doinit is set to 0, this function only updates r->{Lmove, roiidx, inroi}. r->roitype update is not necessary
 * in the latter case, if reference element is found, it also updates r->refeid and r->roisize pointers
 */
void traceroi(ray* r, raytracer* tracer, int roitype, int doinit) {
    int eid = r->eid - 1;
    int* ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);

    if (roitype == 1) { /** edge and node immc - edge also depends on node */
        int i;
        float minratio = 1.f;
        int hitstatus = htNone, firsthit = htNone, firstinout = htNone;

        if (tracer->mesh->edgeroi) { /** if edge roi is defined */
            // edge-based iMMC  - ray-cylinder intersection test
            float distdata[2];
            MMCfloat3 projdata[2];

            for (i = 0; i < 6; i++) { /** loop over each edge in current element */
                if (r->roisize[i] > 0.f) {
                    /** decide if photon is in the roi or not */
                    compute_distances_to_edge(r, tracer, ee, i, distdata, projdata, &hitstatus);

                    /**
                       hitstatus has 4 possible outputs:
                       htInOut: photon path intersects with cylinder, moving from in to out
                       htOutIn: photon path intersects with cylinder, moving from out to in
                       htNoHitIn: both ends of photon path are inside cylinder, no intersection
                       htNoHitOut: both ends of photon path are outside cylinder, no intersection
                       htNone: unexpected, should never happen
                     */
                    if (doinit) {
                        r->inroi |= (hitstatus == htInOut || hitstatus == htNoHitIn);
                    } else {
                        if (hitstatus == htInOut || hitstatus == htOutIn) { /** if intersection is found */
                            /** calculate the first intersection distance normalied by path seg length */
                            float lratio = ray_cylinder_intersect(r, i, distdata, projdata, hitstatus);

                            if (lratio < minratio) {
                                minratio = lratio;
                                firsthit = hitstatus;
                                r->roiidx = i;
                            }
                        } else if (hitstatus == htNoHitIn || hitstatus == htNoHitOut) {
                            firstinout = ((firstinout == htNone) ? hitstatus : (hitstatus == htNoHitIn ? hitstatus : firstinout));
                        }
                    }
                }
            }

            if (minratio < 1.f) {
                r->Lmove *= minratio;
            }

            if (!doinit) {
                r->inroi = (firsthit != htNone ? (firsthit == htOutIn) : (firstinout != htNone ? (firstinout == htNoHitIn) : r->inroi ));
                r->inroi = (firsthit == htNone && firstinout == htNone) ? 0 : r->inroi;
            }

            r->roitype = (firsthit == htInOut || firsthit == htOutIn) ? rtEdge : rtNone;
        }

        if (firsthit == htNone && firstinout != htNoHitIn && tracer->mesh->noderoi) {
            // not hit any edgeroi in the current element, then go for node-based iMMC
            float nr;
            MMCfloat3* center;

            minratio = r->Lmove;
            hitstatus = htNone;
            firsthit = htNone;
            firstinout = htNone;

            for (i = 0; i < 4; i++) { // check if hits any node edgeroi
                nr = tracer->mesh->noderoi[ee[i] - 1];

                if (nr > 0.f) {
                    compute_distances_to_node(r, tracer, ee, i, nr, &center, &hitstatus);

                    if (doinit) {
                        r->inroi |= (hitstatus == htInOut || hitstatus == htNoHitIn);
                    } else {
                        if (hitstatus == htInOut || hitstatus == htOutIn) {
                            float lmove = ray_sphere_intersect(r, i, center, nr, hitstatus);

                            if (lmove < minratio) {
                                minratio = lmove;
                                firsthit = hitstatus;
                                r->roiidx = i;
                            }
                        } else if (hitstatus == htNoHitIn || hitstatus == htNoHitOut) {
                            firstinout = ((firstinout == htNone) ? hitstatus : (hitstatus == htNoHitIn ? hitstatus : firstinout));
                        }
                    }
                }
            }

            if (minratio < r->Lmove) {
                r->Lmove = minratio;
            }

            if (!doinit) {
                r->inroi = (firsthit != htNone ? (firsthit == htOutIn) : (firstinout != htNone ? (firstinout == htNoHitIn) : r->inroi ));
                r->inroi = (firsthit == htNone && firstinout == htNone) ? 0 : r->inroi;
            }

            r->roitype = (firsthit == htInOut || firsthit == htOutIn) ? rtNode : rtNone;
        }

    } else if (tracer->mesh->faceroi) {
        int neweid = -1, newbaseid = 0;
        int* newee;
        float lratio = 1.f, minratio = 1.f;
        int hitstatus = htNone, firsthit = htNone, firstinout = htNone, i;

        // test if this is a reference element, indicated by a negative radius
        if (r->roisize[0] < -4.f) {
            neweid = (int)(-r->roisize[0]) - 4;
            r->refeid = neweid;
            newbaseid = (neweid - 1) << 2;
            newee = (int*)(tracer->mesh->elem + (neweid - 1) * tracer->mesh->elemlen);
            r->roisize = (float*)(tracer->mesh->faceroi + (neweid - 1) * 4);
        }

        // test if the ray intersect with a face ROI (slab)
        for (i = 0; i < 4; i++) {
            if (r->roisize[i] > 0.f) {
                if (neweid < 0) {
                    lratio = ray_face_intersect(r, tracer, ee, i, eid << 2, eid, &hitstatus);
                } else {
                    lratio = ray_face_intersect(r, tracer, newee, i, newbaseid, neweid, &hitstatus);
                }

                if (doinit) {
                    r->inroi |= (hitstatus == htInOut || hitstatus == htNoHitIn);
                } else {
                    if (hitstatus == htInOut || hitstatus == htOutIn) { /** if intersection is found */
                        if (lratio < minratio) {
                            minratio = lratio;
                            firsthit = hitstatus;
                            r->roiidx = i;
                        }
                    } else if (hitstatus == htNoHitIn || hitstatus == htNoHitOut) {
                        firstinout = ((firstinout == htNone) ? hitstatus : (hitstatus == htNoHitIn ? hitstatus : firstinout));
                    }
                }
            }
        }

        if (minratio < 1.f) {
            r->Lmove *= minratio;
        }

        if (!doinit) {
            r->inroi = (firsthit != htNone ? (firsthit == htOutIn) : (firstinout != htNone ? (firstinout == htNoHitIn) : r->inroi ));
            r->inroi = (firsthit == htNone && firstinout == htNone) ? 0 : r->inroi;
        }

        r->roitype = (firsthit == htInOut || firsthit == htOutIn) ? rtFace : rtNone;
    }
}

/**
 * \brief Branch-less Badouel-based SSE4 ray-tracer to advance photon by one step
 *
 * this function uses Branch-less Badouel-based SSE4 ray-triangle intersection
 * tests to advance photon by one step, see Fang2012. Both Badouel and
 * Branch-less Badouel algorithms do not calculate the Barycentric coordinates
 * and can only store energy loss using 0-th order basis function. This function
 * is the fastest among the 4 ray-tracers.
 *
 * \param[in,out] r: the current ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

float branchless_badouel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {

    MMCfloat3 bary = {1e10f, 0.f, 0.f, 0.f};
    float Lp0 = 0.f, rc, currweight, dlen, ww, totalloss = 0.f;
    int tshift, faceidx = -1, baseid, eid;
    __m128 O, T, S;
    __m128i P;

    if (tracer->mesh == NULL || tracer->n == NULL || r->eid <= 0 || r->eid > tracer->mesh->ne) {
        return -1;
    }

    r->p0.w = 1.f;
    r->vec.w = 0.f;
    eid = r->eid - 1;
    baseid = eid << 2;

    r->pout.x = MMC_UNDEFINED;
    r->faceid = -1;
    r->isend = 0;
    r->roitype = rtNone;
    r->refeid = -1;
    r->roiidx = -1;

    const __m128 Nx = _mm_load_ps(&(tracer->n[baseid].x));
    const __m128 Ny = _mm_load_ps(&(tracer->n[baseid + 1].x));
    const __m128 Nz = _mm_load_ps(&(tracer->n[baseid + 2].x));
    const __m128 dd = _mm_load_ps(&(tracer->n[baseid + 3].x));

    O = _mm_set1_ps(r->p0.x);
    T = _mm_mul_ps(Nx, O);
    O = _mm_set1_ps(r->p0.y);
    T = _mm_add_ps(T, _mm_mul_ps(Ny, O));
    O = _mm_set1_ps(r->p0.z);
    T = _mm_add_ps(T, _mm_mul_ps(Nz, O));
    T = _mm_sub_ps(dd, T);

    O = _mm_set1_ps(r->vec.x);
    S = _mm_mul_ps(Nx, O);
    O = _mm_set1_ps(r->vec.y);
    S = _mm_add_ps(S, _mm_mul_ps(Ny, O));
    O = _mm_set1_ps(r->vec.z);
    S = _mm_add_ps(S, _mm_mul_ps(Nz, O));
    T = _mm_div_ps(T, S);

    O = _mm_cmpgt_ps(S, _mm_set1_ps(0.f));
    T = _mm_add_ps(_mm_andnot_ps(O, _mm_set1_ps(1e10f)), _mm_and_ps(O, T));
    S = _mm_movehl_ps(T, T);
    O = _mm_min_ps(T, S);
    S = _mm_shuffle_ps(O, O, _MM_SHUFFLE(1, 1, 1, 1));
    O = _mm_min_ss(O, S);

    _mm_store_ss(&(bary.x), O);
    faceidx = maskmap[_mm_movemask_ps(_mm_cmpeq_ps(T, _mm_set1_ps(bary.x)))];
    r->faceid = faceorder[faceidx];

    if (r->faceid >= 0 && bary.x >= 0) {
        medium* prop;
        int* enb, *ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);
        float mus;

        if (cfg->implicit && r->inroi) {
            prop = tracer->mesh->med + tracer->mesh->prop;
        } else {
            prop = tracer->mesh->med + (tracer->mesh->type[eid]);
        }

        rc = prop->n * R_C0;
        currweight = r->weight;
        mus = (cfg->mcmethod == mmMCX) ? prop->mus : (prop->mua + prop->mus);

        enb = (int*)(tracer->mesh->facenb + eid * tracer->mesh->elemlen);
        r->nexteid = enb[r->faceid]; // if I use nexteid-1, the speed got slower, strange!

        if (r->nexteid > 0) {
            _mm_prefetch((char*) & (tracer->n[(r->nexteid - 1) << 2].x), _MM_HINT_T0);
        }

        dlen = (mus <= EPS) ? R_MIN_MUS : r->slen / mus;
        Lp0 = bary.x;
        r->isend = (Lp0 > dlen);
        r->Lmove = ((r->isend) ? dlen : Lp0);

        // implicit MMC - test if ray intersects with edge/face/node ROI boundaries
        if (cfg->implicit) {
            traceroi(r, tracer, cfg->implicit, 0);
        }

        O = _mm_load_ps(&(r->vec.x));
        S = _mm_load_ps(&(r->p0.x));
        T = _mm_set1_ps(bary.x);
        T = _mm_mul_ps(O, T);
        T = _mm_add_ps(T, S);
        _mm_store_ps(&(r->pout.x), T);

        if ((int)((r->photontimer + r->Lmove * rc - cfg->tstart)*visit->rtstep) > cfg->maxgate - 1) { /*exit time window*/
            r->faceid = -2;
            r->pout.x = MMC_UNDEFINED;
            r->Lmove = (cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f;
        }

        if (cfg->mcmethod == mmMCX) {
            totalloss = expf(-prop->mua * r->Lmove);
            r->weight *= totalloss;
        }

        totalloss = 1.f - totalloss;

        if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otJacobian) {
            currweight = expf(-DELTA_MUA * r->Lmove);
            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        } else if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWL) {
            currweight = r->Lmove;
            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        }

        r->slen -= r->Lmove * mus;

        if (cfg->seed == SEED_FROM_FILE && cfg->outputtype == otWP) {
            if (r->slen0 < EPS) {
                currweight = 1;
            } else {
                currweight = r->Lmove * mus / r->slen0;
            }

            currweight *= cfg->replayweight[r->photonid];
            currweight += r->weight;
        }

        if (bary.x >= 0.f) {
            int framelen = (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);

            if (cfg->method == rtBLBadouelGrid) {
                framelen = cfg->crop0.z;
            }

            ww = currweight - r->weight;
            r->photontimer += r->Lmove * rc;

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * framelen;
            } else {
                tshift = MIN( ((int)((r->photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * framelen;
            }

            if (cfg->debuglevel & dlAccum) MMC_FPRINTF(cfg->flog, "A %f %f %f %e %d %e\n",
                        r->p0.x, r->p0.y, r->p0.z, bary.x, eid + 1, dlen);

            if (prop->mua > 0.f) {
                r->Eabsorb += ww;

                if (cfg->outputtype != otEnergy && cfg->outputtype != otWP) {
                    ww /= prop->mua;
                }
            }

            T = _mm_set1_ps(r->Lmove);
            T = _mm_add_ps(S, _mm_mul_ps(O, T));
            _mm_store_ps(&(r->p0.x), T);

            if (cfg->mcmethod == mmMCX) {
                if (!cfg->basisorder) {
                    if (cfg->method == rtBLBadouel) {
                        unsigned int newidx = eid + tshift;
                        r->oldidx = (r->oldidx == 0xFFFFFFFF) ? newidx : r->oldidx;

                        if (newidx != r->oldidx) {
                            if (cfg->srctype != stPattern) {
                                if (cfg->isatomic)
                                    #pragma omp atomic
                                    tracer->mesh->weight[r->oldidx] += r->oldweight;
                                else {
                                    tracer->mesh->weight[r->oldidx] += r->oldweight;
                                }
                            } else {
                                int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                                int pidx; // pattern index

                                for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                    if (cfg->isatomic)
                                        #pragma omp atomic
                                        tracer->mesh->weight[r->oldidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[pidx * psize + r->posidx];
                                    else {
                                        tracer->mesh->weight[r->oldidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[pidx * psize + r->posidx];
                                    }
                                }
                            }

                            r->oldidx = newidx;
                            r->oldweight = ww;
                        } else {
                            r->oldweight += ww;
                        }

                        if (r->faceid == -2 || !r->isend) {
                            if (cfg->srctype != stPattern) {
                                if (cfg->isatomic)
                                    #pragma omp atomic
                                    tracer->mesh->weight[newidx] += r->oldweight;
                                else {
                                    tracer->mesh->weight[newidx] += r->oldweight;
                                }
                            } else {
                                int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                                int pidx; // pattern index

                                for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                    if (cfg->isatomic)
                                        #pragma omp atomic
                                        tracer->mesh->weight[newidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[pidx * psize + r->posidx];
                                    else {
                                        tracer->mesh->weight[newidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[pidx * psize + r->posidx];
                                    }
                                }
                            }

                            r->oldweight = 0.f;
                        }
                    } else {
                        float dstep, segloss, w0;
                        int4 idx __attribute__ ((aligned(16)));
                        int i, seg = (int)(r->Lmove / cfg->steps.x) + 1;
                        seg = (seg << 1);
                        dstep = r->Lmove / seg;
                        segloss = expf(-prop->mua * dstep);
                        T =  _mm_mul_ps(O, _mm_set1_ps(dstep)); /*step*/
                        O =  _mm_sub_ps(S, _mm_load_ps(&(tracer->mesh->nmin.x)));
                        S =  _mm_add_ps(O, _mm_mul_ps(T, _mm_set1_ps(0.5f))); /*starting point*/
                        dstep = 1.f / cfg->steps.x;
                        totalloss = (totalloss == 0.f) ? 0.f : (1.f - segloss) / totalloss;
                        w0 = ww;

                        for (i = 0; i < seg; i++) {
                            P = _mm_cvttps_epi32(_mm_mul_ps(S, _mm_set1_ps(dstep)));
                            _mm_store_si128((__m128i*) & (idx.x), P);
                            unsigned int newidx = idx.z * cfg->crop0.y + idx.y * cfg->crop0.x + idx.x + tshift;
                            r->oldidx = (r->oldidx == 0xFFFFFFFF) ? newidx : r->oldidx;

                            if (newidx != r->oldidx) {
                                if (cfg->isatomic)
                                    #pragma omp atomic
                                    tracer->mesh->weight[r->oldidx] += r->oldweight;
                                else {
                                    tracer->mesh->weight[r->oldidx] += r->oldweight;
                                }

                                r->oldidx = newidx;
                                r->oldweight = w0 * totalloss;
                            } else {
                                r->oldweight += w0 * totalloss;
                            }

                            if (r->faceid == -2 || !r->isend) {
                                if (cfg->isatomic)
                                    #pragma omp atomic
                                    tracer->mesh->weight[newidx] += r->oldweight;
                                else {
                                    tracer->mesh->weight[newidx] += r->oldweight;
                                }

                                r->oldweight = 0.f;
                            }

                            w0 *= segloss;
                            S = _mm_add_ps(S, T);
                        }
                    }
                } else {
                    int i;
                    ww *= 1.f / 3.f;

                    if (cfg->srctype != stPattern) {
                        if (cfg->isatomic)
                            for (i = 0; i < 3; i++)
                                #pragma omp atomic
                                tracer->mesh->weight[ee[out[faceidx][i]] - 1 + tshift] += ww;
                        else
                            for (i = 0; i < 3; i++) {
                                tracer->mesh->weight[ee[out[faceidx][i]] - 1 + tshift] += ww;
                            }
                    } else {
                        int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w; // total number of pixels in each pattern
                        int pidx; // pattern index

                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                            if (cfg->isatomic)
                                for (i = 0; i < 3; i++)
                                    #pragma omp atomic
                                    tracer->mesh->weight[(ee[out[faceidx][i]] - 1 + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx];
                            else
                                for (i = 0; i < 3; i++) {
                                    tracer->mesh->weight[(ee[out[faceidx][i]] - 1 + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[pidx * psize + r->posidx];
                                }
                        }
                    }
                }
            }
        }
    }

    visit->raytet++;

    if (tracer->mesh->type[eid] == 0) {
        visit->raytet0++;
    }

    r->p0.w = 0.f;

    if (r->faceid == -2) {
        return 0.f;
    }

    return r->slen;
}
#else

/**
 * Havel, Badouel and Branch-less Badouel-based SSE4 ray-tracer require to compile
 * with SSE4 only. Give an error if not compiled with SSE.
 */


float havel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {
    MMC_ERROR(-6, "wrong option, please recompile with SSE4 enabled");
    return MMC_UNDEFINED;
}
float badouel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {
    MMC_ERROR(-6, "wrong option, please recompile with SSE4 enabled");
    return MMC_UNDEFINED;
}
float branchless_badouel_raytet(ray* r, raytracer* tracer, mcconfig* cfg, visitor* visit) {
    MMC_ERROR(-6, "wrong option, please recompile with SSE4 enabled");
    return MMC_UNDEFINED;
}
#endif

/**
 * @brief The core Monte Carlo function simulating a single photon (!!!Important!!!)
 *
 * This is the core Monte Carlo simulation function. It simulates the life-time
 * of a single photon packet, from launching to termination.
 *
 * \param[in] id: the linear index of the current photon, starting from 0.
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] mesh: the mesh data structure
 * \param[in,out] ran: the random number generator states
 * \param[in,out] ran0: the additional random number generator states
 * \param[in,out] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

void onephoton(size_t id, raytracer* tracer, tetmesh* mesh, mcconfig* cfg,
               RandType* ran, RandType* ran0, visitor* visit) {

    int oldeid, fixcount = 0, exitdet = 0;
    int* enb;
    float mom;
    float kahany, kahant;
    ray r = {cfg->srcpos, {cfg->srcdir.x, cfg->srcdir.y, cfg->srcdir.z}, {MMC_UNDEFINED, 0.f, 0.f}, cfg->bary0, cfg->e0, cfg->dim.y - 1, 0, 0, 1.f, 0.f, 0.f, 0.f, 0.f, 0., 0, NULL, NULL, cfg->srcdir.w, 0, 0xFFFFFFFF, 0.0, NULL, 0, 0, 0, 0};

    float (*engines[5])(ray * r, raytracer * tracer, mcconfig * cfg, visitor * visit) =
    {plucker_raytet, havel_raytet, badouel_raytet, branchless_badouel_raytet, branchless_badouel_raytet};
    float (*tracercore)(ray * r, raytracer * tracer, mcconfig * cfg, visitor * visit);

    r.partialpath = (float*)calloc(visit->reclen - 1, sizeof(float));
    r.photonid = id;

    if (cfg->issavedet && cfg->issaveseed) {
        r.photonseed = (void*)calloc(1, (sizeof(RandType) * RAND_BUF_LEN));
        memcpy(r.photonseed, (void*)ran, (sizeof(RandType)*RAND_BUF_LEN));
    }

    tracercore = engines[0];

    if (cfg->method >= rtPlucker && cfg->method <= rtBLBadouelGrid) {
        tracercore = engines[(int)(cfg->method)];
    } else {
        MMC_ERROR(-6, "specified ray-tracing algorithm is not defined");
    }

    /*initialize the photon parameters*/
    launchphoton(cfg, &r, mesh, ran, ran0);
    r.partialpath[visit->reclen - 2] = r.weight; /*last record in partialpath is the initial photon weight*/

    /*use Kahan summation to accumulate weight, otherwise, counter stops at 16777216*/
    /*http://stackoverflow.com/questions/2148149/how-to-sum-a-large-number-of-float-number*/
    int pidx;

    if (cfg->srctype != stPattern) {
        if (cfg->seed == SEED_FROM_FILE && (cfg->outputtype == otWL || cfg->outputtype == otWP)) {
            kahany = cfg->replayweight[r.photonid] - visit->kahanc0[0];    /* when replay mode, accumulate detected photon weight */
        } else {
            kahany = r.weight - visit->kahanc0[0];
        }

        kahant = visit->launchweight[0] + kahany;
        visit->kahanc0[0] = (kahant - visit->launchweight[0]) - kahany;
        visit->launchweight[0] = kahant;
    } else {
        int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w;
        #pragma omp critical
        {
            for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                if (cfg->seed == SEED_FROM_FILE && (cfg->outputtype == otWL || cfg->outputtype == otWP)) {
                    kahany = cfg->replayweight[r.photonid] - visit->kahanc0[pidx];    /* when replay mode, accumulate detected photon weight */
                } else {
                    kahany = r.weight * cfg->srcpattern[pidx * psize + r.posidx] - visit->kahanc0[pidx];
                }

                kahant = visit->launchweight[pidx] + kahany;
                visit->kahanc0[pidx] = (kahant - visit->launchweight[pidx]) - kahany;
                visit->launchweight[pidx] = kahant;
            }
        }
    }

#ifdef MMC_USE_SSE
    const float int_coef_arr[4] = { -1.f, -1.f, -1.f, 1.f };
    int_coef = _mm_load_ps(int_coef_arr);

    if (cfg->implicit) {
        updateroi(cfg->implicit, &r, tracer->mesh);
        traceroi(&r, tracer, cfg->implicit, 1);
    }

#endif

    while (1) { /*propagate a photon until exit*/
        if (cfg->implicit) {
            updateroi(cfg->implicit, &r, tracer->mesh);
        }

        r.slen = (*tracercore)(&r, tracer, cfg, visit);

        if (r.pout.x == MMC_UNDEFINED) {
            if (r.faceid == -2) {
                break;    /*reaches the time limit*/
            }

            if (fixcount++ < MAX_TRIAL) {
                fixphoton(&r.p0, mesh->node, (int*)(mesh->elem + (r.eid - 1)*mesh->elemlen));
                continue;
            }

            r.eid = ID_UNDEFINED;
            r.faceid = -1;
        }

        if (cfg->issavedet && r.Lmove > 0.f && mesh->type[r.eid - 1] > 0 && r.faceid >= 0) {
            r.partialpath[mesh->prop - 1 + mesh->type[r.eid - 1]] += r.Lmove;    /*second medianum block is the partial path*/
        }

        if (cfg->implicit && cfg->isreflect && r.roitype && r.roiidx >= 0 && (mesh->med[cfg->his.maxmedia].n != mesh->med[mesh->type[r.eid - 1]].n)) {
            reflectrayroi(cfg, &r.vec, &r.p0, tracer, &r.eid, &r.inroi, ran, r.roitype, r.roiidx, r.refeid);
            vec_mult_add(&r.p0, &r.vec, 1.0f, 10 * EPS, &r.p0);
            continue;
        } else if (cfg->implicit && r.roitype) {
            continue;
        }

        /*move a photon until the end of the current scattering path*/
        while (r.faceid >= 0 && !r.isend && !r.roitype) {
            memcpy((void*)&r.p0, (void*)&r.pout, sizeof(r.p0));

            enb = (int*)(mesh->facenb + (r.eid - 1) * mesh->elemlen);
            oldeid = r.eid;
            r.eid = enb[r.faceid];

            if (cfg->implicit) {
                if (cfg->isreflect && (r.eid <= 0 || mesh->med[mesh->type[r.eid - 1]].n != mesh->med[mesh->type[oldeid - 1]].n )) {
                    if (! (!r.inroi && r.eid <= 0 && ((mesh->med[mesh->type[oldeid - 1]].n == cfg->nout && cfg->isreflect != (int)bcMirror) || cfg->isreflect == (int)bcAbsorbExterior) ) ) {
                        reflectray(cfg, &r.vec, tracer, &oldeid, &r.eid, r.faceid, ran, r.inroi);
                    }
                }
            } else {
                if (cfg->isreflect && (r.eid <= 0 || mesh->med[mesh->type[r.eid - 1]].n != mesh->med[mesh->type[oldeid - 1]].n )) {
                    if (! (r.eid <= 0 && ((mesh->med[mesh->type[oldeid - 1]].n == cfg->nout && cfg->isreflect != (int)bcMirror) || cfg->isreflect == (int)bcAbsorbExterior) ) ) {
                        reflectray(cfg, &r.vec, tracer, &oldeid, &r.eid, r.faceid, ran, r.inroi);
                    }
                }
            }

            if (r.eid <= 0) {
                break;
            }

            /*when a photon enters the domain from the background*/
            if (mesh->type[oldeid - 1] == 0 && mesh->type[r.eid - 1]) {
                if (cfg->debuglevel & dlExit)
                    MMC_FPRINTF(cfg->flog, "e %f %f %f %f %f %f %f %d\n", r.p0.x, r.p0.y, r.p0.z,
                                r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid);

                if (!cfg->voidtime) {
                    r.photontimer = 0.f;
                }
            }

            /*when a photon exits the domain into the background*/
            if (mesh->type[oldeid - 1] && mesh->type[r.eid - 1] == 0) {
                if (cfg->debuglevel & dlExit)
                    MMC_FPRINTF(cfg->flog, "x %f %f %f %f %f %f %f %d\n", r.p0.x, r.p0.y, r.p0.z,
                                r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid);

                if (!cfg->isextdet) {
                    r.eid = 0;
                    break;
                }
            }

            //          if(r.eid!=ID_UNDEFINED && mesh->med[mesh->type[oldeid-1]].n == cfg->nout ) break;
            if (r.pout.x != MMC_UNDEFINED && (cfg->debuglevel & dlMove)) {
                MMC_FPRINTF(cfg->flog, "P %f %f %f %d %zu %e\n", r.pout.x, r.pout.y, r.pout.z, r.eid, id, r.slen);
            }

            if (cfg->implicit) {
                updateroi(cfg->implicit, &r, tracer->mesh);
            }

            r.slen = (*tracercore)(&r, tracer, cfg, visit);

            if (cfg->issavedet && r.Lmove > 0.f && mesh->type[r.eid - 1] > 0) {
                r.partialpath[mesh->prop - 1 + mesh->type[r.eid - 1]] += r.Lmove;
            }

            if (r.faceid == -2) {
                break;
            }

            fixcount = 0;

            while (r.pout.x == MMC_UNDEFINED && fixcount++ < MAX_TRIAL) {
                fixphoton(&r.p0, mesh->node, (int*)(mesh->elem + (r.eid - 1)*mesh->elemlen));
                r.slen = (*tracercore)(&r, tracer, cfg, visit);

                if (cfg->issavedet && r.Lmove > 0.f && mesh->type[r.eid - 1] > 0) {
                    r.partialpath[mesh->prop - 1 + mesh->type[r.eid - 1]] += r.Lmove;
                }
            }

            if (r.pout.x == MMC_UNDEFINED) {
                /*possibily hit an edge or miss*/
                r.eid = ID_UNDEFINED;
                break;
            }
        }

        if (r.eid <= 0 || r.pout.x == MMC_UNDEFINED) {
            if (r.eid != ID_UNDEFINED && (cfg->debuglevel & dlMove)) {
                MMC_FPRINTF(cfg->flog, "B %f %f %f %d %zu %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen);
            }

            if (r.eid != ID_UNDEFINED) {
                if (cfg->debuglevel & dlExit)
                    MMC_FPRINTF(cfg->flog, "E %f %f %f %f %f %f %f %d\n", r.p0.x, r.p0.y, r.p0.z,
                                r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid);

                if (cfg->issavedet && cfg->issaveexit) {                                   /*when issaveexit is set to 1*/
                    memcpy(r.partialpath + (visit->reclen - 2 - 6), &(r.p0.x), sizeof(float) * 3); /*columns 7-5 from the right store the exit positions*/
                    memcpy(r.partialpath + (visit->reclen - 2 - 3), &(r.vec.x), sizeof(float) * 3); /*columns 4-2 from the right store the exit dirs*/
                }

                if (cfg->issaveref && r.eid < 0 && mesh->dref) {
                    int tshift = MIN( ((int)((r.photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * mesh->nf;
                    mesh->dref[((-r.eid) - 1) + tshift] += r.weight;
                }
            } else if (r.faceid == -2 && (cfg->debuglevel & dlMove)) {
                MMC_FPRINTF(cfg->flog, "T %f %f %f %d %zu %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen);
            } else if (r.eid && r.faceid != -2  && cfg->debuglevel & dlEdge) {
                MMC_FPRINTF(cfg->flog, "X %f %f %f %d %zu %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen);
            }

            if (cfg->issavedet && r.eid <= 0) {
                int i;

                if (cfg->detnum == 0 && cfg->isextdet && mesh->type[oldeid - 1] == mesh->prop + 1) {
                    exitdet = oldeid;
                } else
                    for (i = 0; i < cfg->detnum; i++) {
                        if ((cfg->detpos[i].x - r.p0.x) * (cfg->detpos[i].x - r.p0.x) +
                                (cfg->detpos[i].y - r.p0.y) * (cfg->detpos[i].y - r.p0.y) +
                                (cfg->detpos[i].z - r.p0.z) * (cfg->detpos[i].z - r.p0.z) < cfg->detpos[i].w * cfg->detpos[i].w) {
                            exitdet = i + 1;
                            break;
                        }
                    }
            }

            break;  /*photon exits boundary*/
        }

        if (cfg->debuglevel & dlMove) {
            MMC_FPRINTF(cfg->flog, "M %f %f %f %d %zu %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen);
        }

        if (cfg->minenergy > 0.f && r.weight < cfg->minenergy && (cfg->tend - cfg->tstart)*visit->rtstep <= 1.f) { /*Russian Roulette*/
            if (rand_do_roulette(ran)*cfg->roulettesize <= 1.f) {
                r.weight *= cfg->roulettesize;

                if (cfg->debuglevel & dlWeight) {
                    MMC_FPRINTF(cfg->flog, "Russian Roulette bumps r.weight to %f\n", r.weight);
                }
            } else {
                break;
            }
        }

        if (cfg->implicit && cfg->isreflect && r.roitype && r.roiidx >= 0 && (mesh->med[cfg->his.maxmedia].n != mesh->med[mesh->type[r.eid - 1]].n)) {
            reflectrayroi(cfg, &r.vec, &r.p0, tracer, &r.eid, &r.inroi, ran, r.roitype, r.roiidx, r.refeid);
            vec_mult_add(&r.p0, &r.vec, 1.0f, 10 * EPS, &r.p0);
            continue;
        } else if (cfg->implicit && r.roitype) {
            continue;
        }

        mom = 0.f;
        r.slen0 = mc_next_scatter(mesh->med[mesh->type[r.eid - 1]].g, &r.vec, ran, ran0, cfg, &mom);
        r.slen = r.slen0;

        if (cfg->mcmethod != mmMCX) {
            albedoweight(&r, mesh, cfg, visit);
        }

        if (cfg->ismomentum && mesh->type[r.eid - 1] > 0) {              /*when ismomentum is set to 1*/
            r.partialpath[(mesh->prop << 1) - 1 + mesh->type[r.eid - 1]] += mom;    /*the third medianum block stores the momentum transfer*/
        }

        r.partialpath[mesh->type[r.eid - 1] - 1]++;                      /*the first medianum block stores the scattering event counts*/
    }

    if (cfg->issavedet && exitdet > 0) {
        int offset = visit->bufpos * visit->reclen;

        if (visit->bufpos >= visit->detcount) {
            visit->detcount += DET_PHOTON_BUF;
            visit->partialpath = (float*)realloc(visit->partialpath,
                                                 visit->detcount * visit->reclen * sizeof(float));

            if (cfg->issaveseed) {
                visit->photonseed = realloc(visit->photonseed, visit->detcount * (sizeof(RandType) * RAND_BUF_LEN));
            }
        }

        visit->partialpath[offset] = exitdet;
        memcpy(visit->partialpath + offset + 1, r.partialpath, (visit->reclen - 1)*sizeof(float));

        if (cfg->issaveseed) {
            memcpy(visit->photonseed + visit->bufpos * (sizeof(RandType)*RAND_BUF_LEN), r.photonseed, (sizeof(RandType)*RAND_BUF_LEN));
        }

        visit->bufpos++;
    }

    free(r.partialpath);

    if (r.photonseed) {
        free(r.photonseed);
    }

    if (cfg->srctype != stPattern) {
        kahany = r.Eabsorb - visit->kahanc1[0];
        kahant = visit->absorbweight[0] + kahany;
        visit->kahanc1[0] = (kahant - visit->absorbweight[0]) - kahany;
        visit->absorbweight[0] = kahant;
    } else {
        int psize = (int)cfg->srcparam1.w * (int)cfg->srcparam2.w;
        #pragma omp critical
        {
            for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                kahany = r.Eabsorb * cfg->srcpattern[pidx * psize + r.posidx] - visit->kahanc1[pidx];
                kahant = visit->absorbweight[pidx] + kahany;
                visit->kahanc1[pidx] = (kahant - visit->absorbweight[pidx]) - kahany;
                visit->absorbweight[pidx] = kahant;
            }
        }
    }
}

/**
 * @brief Calculate the reflection/transmission of a ray at the ROI surface
 *
 * This function handles the reflection and transmission events at the roi surface wall
 * where the refractive indices mismatch.
 *
 * \param[in] c0: the current direction vector of the ray
 * \param[in] u: the current direction vector of the edge
 * \param[in] ph: hitting position at the edgeroi
 * \param[in] E0: any point on the edge
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] eid: the index of the CURRENT element
 * \param[in,out] ran: the random number generator states
 */

float reflectrayroi(mcconfig* cfg, MMCfloat3* c0, MMCfloat3* ph, raytracer* tracer, int* eid, int* inroi, RandType* ran, int roitype, int roiidx, int refeid) {

    /*to handle refractive index mismatch*/
    MMCfloat3 pnorm = {0.f}, *pn = &pnorm, EH, ut;
    float Icos, Re, Im, Rtotal, tmp0, tmp1, tmp2, n1, n2;

    if (roitype == rtEdge) { // hit edge roi
        MMCfloat3 u, E0;
        E0 = tracer->mesh->node[tracer->mesh->elem[((*eid - 1) << 2) + e2n[roiidx][0]] - 1];
        vec_diff(&E0, tracer->mesh->node + (tracer->mesh->elem[((*eid - 1) << 2) + e2n[roiidx][1]] - 1), &u);
        vec_diff(&E0, ph, &EH);
        tmp0 = vec_dot(&EH, &u);
        vec_mult(&u, tmp0, &ut);
        vec_diff(&ut, &EH, pn);
        tmp0 = 1.f / sqrtf(vec_dot(pn, pn));
        vec_mult(pn, tmp0, pn);
    } else if (roitype == rtNode) { // hit node roi
        vec_diff(tracer->mesh->node + (tracer->mesh->elem[((*eid - 1) << 2) + roiidx] - 1), ph, pn);
        tmp0 = 1.f / sqrtf(vec_dot(pn, pn));
        vec_mult(pn, tmp0, pn);
    } else {         // hit face roi
        int baseid;

        if (refeid < 0) {
            baseid = (*eid - 1) << 2;
        } else {
            baseid = (refeid - 1) << 2;
        }

        pn->x = (&(tracer->n[baseid].x))[roiidx];
        pn->y = (&(tracer->n[baseid].x))[roiidx + 4];
        pn->z = (&(tracer->n[baseid].x))[roiidx + 8];
    }

    /*pn pointing outward*/

    // /*compute the cos of the incidence angle*/
    Icos = fabs(vec_dot(c0, pn));

    if (*inroi) { // out->in
        n1 = tracer->mesh->med[tracer->mesh->type[*eid - 1]].n;
        n2 = tracer->mesh->med[cfg->his.maxmedia].n;

        if (roitype == rtEdge || roitype == rtNode) {
            vec_mult(pn, -1.f, pn);
        }
    } else {     // in->out
        n1 = tracer->mesh->med[cfg->his.maxmedia].n;
        n2 = tracer->mesh->med[tracer->mesh->type[*eid - 1]].n;

        if (roitype == rtFace) {
            vec_mult(pn, -1.f, pn);
        }
    }

    tmp0 = n1 * n1;
    tmp1 = n2 * n2;
    tmp2 = 1.f - tmp0 / tmp1 * (1.f - Icos * Icos); /*1-[n1/n2*sin(si)]^2 = cos(ti)^2*/

    if (tmp2 > 0.f) { /*if no total internal reflection*/
        Re = tmp0 * Icos * Icos + tmp1 * tmp2; /*transmission angle*/
        tmp2 = sqrtf(tmp2); /*to save one sqrt*/
        Im = 2.f * n1 * n2 * Icos * tmp2;
        Rtotal = (Re - Im) / (Re + Im); /*Rp*/
        Re = tmp1 * Icos * Icos + tmp0 * tmp2 * tmp2;
        Rtotal = (Rtotal + (Re - Im) / (Re + Im)) * 0.5f; /*(Rp+Rs)/2*/

        // if(*oldeid==*eid) return Rtotal; /*initial specular reflection*/
        if (rand_next_reflect(ran) <= Rtotal) { /*do reflection*/
            vec_mult_add(pn, c0, -2.f * Icos, 1.f, c0);
            *inroi = !(*inroi);
            //if(cfg->debuglevel&dlReflect) MMC_FPRINTF(cfg->flog,"R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal);
        } else if (cfg->isspecular == 2 && *eid == 0) {
            // if do transmission, but next neighbor is 0, terminate
        } else {                             /*do transmission*/
            vec_mult_add(pn, c0, -Icos, 1.f, c0);
            vec_mult_add(pn, c0, tmp2, n1 / n2, c0);
        }
    } else { /*total internal reflection*/
        vec_mult_add(pn, c0, -2.f * Icos, 1.f, c0);
        *inroi = !(*inroi);
    }

    tmp0 = 1.f / sqrtf(vec_dot(c0, c0));
    vec_mult(c0, tmp0, c0);
    return 1.f;
}

/**
 * @brief Calculate the reflection/transmission of a ray at an interface
 *
 * This function handles the reflection and transmission events at an interface
 * where the refractive indices mismatch.
 *
 * \param[in,out] cfg: simulation configuration structure
 * \param[in] c0: the current direction vector of the ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] oldeid: the index of the element the photon moves away from
 * \param[in] eid: the index of the element the photon about to move into
 * \param[in] faceid: index of the face through which the photon reflects/transmits
 * \param[in,out] ran: the random number generator states
 */

float reflectray(mcconfig* cfg, MMCfloat3* c0, raytracer* tracer, int* oldeid, int* eid, int faceid, RandType* ran, int inroi) {
    /*to handle refractive index mismatch*/
    MMCfloat3 pnorm = {0.f}, *pn = &pnorm;
    float Icos, Re, Im, Rtotal, tmp0, tmp1, tmp2, n1, n2;
    int offs = (*oldeid - 1) << 2;

    faceid = ifaceorder[faceid];

    /*calculate the normal direction of the intersecting triangle*/
    if (cfg->method == rtPlucker) { //Plucker ray-tracing
        pn = tracer->n + (offs) + faceid;
    } else if (cfg->method < rtBLBadouel) {
        pn = tracer->m + (offs + faceid) * 3;
    } else if (cfg->method == rtBLBadouel || cfg->method == rtBLBadouelGrid) {
        pnorm.x = (&(tracer->n[offs].x))[faceid];
        pnorm.y = (&(tracer->n[offs].x))[faceid + 4];
        pnorm.z = (&(tracer->n[offs].x))[faceid + 8];
    }

    /*pn pointing outward*/

    /*compute the cos of the incidence angle*/
    Icos = fabs(vec_dot(c0, pn));

    if (cfg->implicit && inroi != 0) {
        n1 = tracer->mesh->med[cfg->his.maxmedia].n;
        n2 = n1;
    } else {
        n1 = (*oldeid != *eid) ? tracer->mesh->med[tracer->mesh->type[*oldeid - 1]].n : cfg->nout;
        n2 = (*eid > 0) ? tracer->mesh->med[tracer->mesh->type[*eid - 1]].n : cfg->nout;
    }

    tmp0 = n1 * n1;
    tmp1 = n2 * n2;
    tmp2 = 1.f - tmp0 / tmp1 * (1.f - Icos * Icos); /*1-[n1/n2*sin(si)]^2 = cos(ti)^2*/

    if (tmp2 > 0.f && !(*eid <= 0 && cfg->isreflect == bcMirror)) { /*if no total internal reflection*/
        Re = tmp0 * Icos * Icos + tmp1 * tmp2; /*transmission angle*/
        tmp2 = sqrtf(tmp2); /*to save one sqrt*/
        Im = 2.f * n1 * n2 * Icos * tmp2;
        Rtotal = (Re - Im) / (Re + Im); /*Rp*/
        Re = tmp1 * Icos * Icos + tmp0 * tmp2 * tmp2;
        Rtotal = (Rtotal + (Re - Im) / (Re + Im)) * 0.5f; /*(Rp+Rs)/2*/

        if (*oldeid == *eid) {
            return Rtotal;    /*initial specular reflection*/
        }

        if (rand_next_reflect(ran) <= Rtotal) { /*do reflection*/
            vec_mult_add(pn, c0, -2.f * Icos, 1.f, c0);
            //if(cfg->debuglevel&dlReflect) MMC_FPRINTF(cfg->flog,"R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal);
            *eid = *oldeid; /*stay with the current element*/
        } else if (cfg->isspecular == 2 && *eid == 0) {
            // if do transmission, but next neighbor is 0, terminate
        } else {                             /*do transmission*/
            vec_mult_add(pn, c0, -Icos, 1.f, c0);
            vec_mult_add(pn, c0, tmp2, n1 / n2, c0);
            //if(cfg->debuglevel&dlReflect) MMC_FPRINTF(cfg->flog,"Z %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f-Rtotal);
        }
    } else { /*total internal reflection*/
        vec_mult_add(pn, c0, -2.f * Icos, 1.f, c0);
        *eid = *oldeid;
        //if(cfg->debuglevel&dlReflect) MMC_FPRINTF(cfg->flog,"V %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f);
    }

    tmp0 = 1.f / sqrtf(vec_dot(c0, c0));
    vec_mult(c0, tmp0, c0);
    return 1.f;
}

/**
 * @brief Launch a new photon
 *
 * This function launch a new photon using one of the dozen supported source forms.
 *
 * \param[in,out] cfg: simulation configuration structure
 * \param[in,out] r: the current ray
 * \param[in] mesh: the mesh data structure
 * \param[in,out] ran: the random number generator states
 * \param[in,out] ran0: the additional random number generator states
 */

void launchphoton(mcconfig* cfg, ray* r, tetmesh* mesh, RandType* ran, RandType* ran0) {
    int canfocus = 1;
    MMCfloat3 origin = {r->p0.x, r->p0.y, r->p0.z};

    r->slen = rand_next_scatlen(ran);
    r->inroi = 0;

    if (cfg->srctype == stPencil) { // pencil beam, use the old workflow, except when eid is not given
        if (r->eid > 0) {
            return;
        }
    } else if (cfg->srctype == stPlanar || cfg->srctype == stPattern || cfg->srctype == stFourier) {
        float rx = rand_uniform01(ran);
        float ry = rand_uniform01(ran);
        r->p0.x = cfg->srcpos.x + rx * cfg->srcparam1.x + ry * cfg->srcparam2.x;
        r->p0.y = cfg->srcpos.y + rx * cfg->srcparam1.y + ry * cfg->srcparam2.y;
        r->p0.z = cfg->srcpos.z + rx * cfg->srcparam1.z + ry * cfg->srcparam2.z;
        r->weight = 1.f;

        if (cfg->srctype == stPattern) {
            int xsize = (int)cfg->srcparam1.w;
            int ysize = (int)cfg->srcparam2.w;
            r->posidx = MIN((int)(ry * ysize), ysize - 1) * xsize + MIN((int)(rx * xsize), xsize - 1);

            if (cfg->seed == SEED_FROM_FILE && (cfg->outputtype == otWL || cfg->outputtype == otWP)) { // replay mode currently doesn't support multiple source patterns
                r->weight = cfg->srcpattern[MIN( (int)(ry * cfg->srcparam2.w), (int)cfg->srcparam2.w - 1 ) * (int)(cfg->srcparam1.w) + MIN( (int)(rx * cfg->srcparam1.w), (int)cfg->srcparam1.w - 1 )];
                cfg->replayweight[r->photonid] *= r->weight;
            }
        } else if (cfg->srctype == stFourier) {
            r->weight = (cosf((floorf(cfg->srcparam1.w) * rx + floorf(cfg->srcparam2.w) * ry + cfg->srcparam1.w - floorf(cfg->srcparam1.w)) * TWO_PI) * (1.f - cfg->srcparam2.w + floorf(cfg->srcparam2.w)) + 1.f) * 0.5f;
        }

        origin.x += (cfg->srcparam1.x + cfg->srcparam2.x) * 0.5f;
        origin.y += (cfg->srcparam1.y + cfg->srcparam2.y) * 0.5f;
        origin.z += (cfg->srcparam1.z + cfg->srcparam2.z) * 0.5f;
    } else if (cfg->srctype == stFourierX || cfg->srctype == stFourier2D) {
        float rx = rand_uniform01(ran);
        float ry = rand_uniform01(ran);
        float4 v2 = cfg->srcparam1;
        v2.w *= 1.f / (sqrtf(cfg->srcparam1.x * cfg->srcparam1.x + cfg->srcparam1.y * cfg->srcparam1.y + cfg->srcparam1.z * cfg->srcparam1.z));
        v2.x = v2.w * (cfg->srcdir.y * cfg->srcparam1.z - cfg->srcdir.z * cfg->srcparam1.y);
        v2.y = v2.w * (cfg->srcdir.z * cfg->srcparam1.x - cfg->srcdir.x * cfg->srcparam1.z);
        v2.z = v2.w * (cfg->srcdir.x * cfg->srcparam1.y - cfg->srcdir.y * cfg->srcparam1.x);
        r->p0.x = cfg->srcpos.x + rx * cfg->srcparam1.x + ry * v2.x;
        r->p0.y = cfg->srcpos.y + rx * cfg->srcparam1.y + ry * v2.y;
        r->p0.z = cfg->srcpos.z + rx * cfg->srcparam1.z + ry * v2.z;

        if (cfg->srctype == stFourier2D) {
            r->weight = (sinf((cfg->srcparam2.x * rx + cfg->srcparam2.z) * TWO_PI) * sinf((cfg->srcparam2.y * ry + cfg->srcparam2.w) * TWO_PI) + 1.f) * 0.5f;    //between 0 and 1
        } else {
            r->weight = (cosf((cfg->srcparam2.x * rx + cfg->srcparam2.y * ry + cfg->srcparam2.z) * TWO_PI) * (1.f - cfg->srcparam2.w) + 1.f) * 0.5f;    //between 0 and 1
        }

        origin.x += (cfg->srcparam1.x + v2.x) * 0.5f;
        origin.y += (cfg->srcparam1.y + v2.y) * 0.5f;
        origin.z += (cfg->srcparam1.z + v2.z) * 0.5f;
    } else if (cfg->srctype == stDisk || cfg->srctype == stGaussian) { // uniform disk and Gaussian beam
        float sphi, cphi;
        float phi = TWO_PI * rand_uniform01(ran);
        sphi = sinf(phi);
        cphi = cosf(phi);
        float r0;

        if (cfg->srctype == stDisk) {
            r0 = sqrtf(rand_uniform01(ran)) * cfg->srcparam1.x;
        } else if (fabs(r->focus) < 1e-5f || fabs(cfg->srcparam1.y) < 1e-5f) {
            r0 = sqrtf(-log(rand_uniform01(ran))) * cfg->srcparam1.x;
        } else {
            float z0 = cfg->srcparam1.x * cfg->srcparam1.x * M_PI / cfg->srcparam1.y; //Rayleigh range
            r0 = sqrtf(-log(rand_uniform01(ran)) * (1.f + (r->focus * r->focus / (z0 * z0)))) * cfg->srcparam1.x;
        }

        if (cfg->srcdir.z > -1.f + EPS && cfg->srcdir.z < 1.f - EPS) {
            float tmp0 = 1.f - cfg->srcdir.z * cfg->srcdir.z;
            float tmp1 = r0 / sqrtf(tmp0);
            r->p0.x = cfg->srcpos.x + tmp1 * (cfg->srcdir.x * cfg->srcdir.z * cphi - cfg->srcdir.y * sphi);
            r->p0.y = cfg->srcpos.y + tmp1 * (cfg->srcdir.y * cfg->srcdir.z * cphi + cfg->srcdir.x * sphi);
            r->p0.z = cfg->srcpos.z - tmp1 * tmp0 * cphi;
        } else {
            r->p0.x += r0 * cphi;
            r->p0.y += r0 * sphi;
        }
    } else if (cfg->srctype == stCone || cfg->srctype == stIsotropic || cfg->srctype == stArcSin) {
        float ang, stheta, ctheta, sphi, cphi;
        ang = TWO_PI * rand_uniform01(ran); //next arimuth angle
        sphi = sinf(ang);
        cphi = cosf(ang);

        if (cfg->srctype == stCone) { // a solid-angle section of a uniform sphere
            do {
                ang = (cfg->srcparam1.y > 0) ? TWO_PI * rand_uniform01(ran) : acosf(2.f * rand_uniform01(ran) - 1.f); //sine distribution
            } while (ang > cfg->srcparam1.x);
        } else {
            if (cfg->srctype == stIsotropic) { // uniform sphere
                ang = acosf(2.f * rand_uniform01(ran) - 1.f);    //sine distribution
            } else {
                ang = M_PI * rand_uniform01(ran);    //uniform distribution in zenith angle, arcsine
            }
        }

        stheta = sinf(ang);
        ctheta = cosf(ang);
        r->vec.x = stheta * cphi;
        r->vec.y = stheta * sphi;
        r->vec.z = ctheta;
        canfocus = 0;

        if (r->eid > 0) {
            return;
        }
    } else if (cfg->srctype == stZGaussian) {
        float ang, stheta, ctheta, sphi, cphi;
        ang = TWO_PI * rand_uniform01(ran); //next arimuth angle
        sphi = sinf(ang);
        cphi = cosf(ang);
        ang = sqrtf(-2.f * log(rand_uniform01(ran))) * (1.f - 2.f * rand_uniform01(ran0)) * cfg->srcparam1.x;
        stheta = sinf(ang);
        ctheta = cosf(ang);
        r->vec.x = stheta * cphi;
        r->vec.y = stheta * sphi;
        r->vec.z = ctheta;
        canfocus = 0;
    } else if (cfg->srctype == stLine || cfg->srctype == stSlit) {
        float t = rand_uniform01(ran);
        r->p0.x += t * cfg->srcparam1.x;
        r->p0.y += t * cfg->srcparam1.y;
        r->p0.z += t * cfg->srcparam1.z;

        if (cfg->srctype == stLine) {
            float s, p;
            t = 1.f - 2.f * rand_uniform01(ran);
            s = 1.f - 2.f * rand_uniform01(ran);
            p = sqrtf(1.f - r->vec.x * r->vec.x - r->vec.y * r->vec.y) * (rand_uniform01(ran) > 0.5f ? 1.f : -1.f);
            MMCfloat3 vv;
            vv.x = r->vec.y * p - r->vec.z * s;
            vv.y = r->vec.z * t - r->vec.x * p;
            vv.z = r->vec.x * s - r->vec.y * t;
            r->vec = vv;
            //*((MMCfloat3*)&(r->vec))=(MMCfloat3)(r->vec.y*p-r->vec.z*s,r->vec.z*t-r->vec.x*p,r->vec.x*s-r->vec.y*t);
        }

        origin.x += (cfg->srcparam1.x) * 0.5f;
        origin.y += (cfg->srcparam1.y) * 0.5f;
        origin.z += (cfg->srcparam1.z) * 0.5f;
        canfocus = (cfg->srctype == stSlit);
    }

    if (canfocus && r->focus != 0.f) { // if beam focus is set, determine the incident angle
        float Rn2;
        origin.x += r->focus * r->vec.x;
        origin.y += r->focus * r->vec.y;
        origin.z += r->focus * r->vec.z;

        if (r->focus < 0.f) { // diverging beam
            r->vec.x = r->p0.x - origin.x;
            r->vec.y = r->p0.y - origin.y;
            r->vec.z = r->p0.z - origin.z;
        } else {            // converging beam
            r->vec.x = origin.x - r->p0.x;
            r->vec.y = origin.y - r->p0.y;
            r->vec.z = origin.z - r->p0.z;
        }

        Rn2 = 1.f / sqrtf(vec_dot(&r->vec, &r->vec)); // normalize
        vec_mult(&r->vec, Rn2, &r->vec);
    }

    vec_mult_add(&(r->p0), &(r->vec), 1.f, EPS, &(r->p0));

    /*Caluclate intial element id and bary-centric coordinates for area sources - position changes everytime*/
    MMCfloat3 vecS = {0.f}, *nodes = mesh->node, vecAB, vecAC, vecN;
    int is, i, ea, eb, ec;
    float bary[4] = {0.f};

    for (is = 0; is < mesh->srcelemlen; is++) {
        int include = 1;
        int* elems = (int*)(mesh->elem + (mesh->srcelem[is] - 1) * mesh->elemlen);

        for (i = 0; i < 4; i++) {
            ea = elems[out[i][0]] - 1;
            eb = elems[out[i][1]] - 1;
            ec = elems[out[i][2]] - 1;
            vec_diff(&nodes[ea], &nodes[eb], &vecAB); //repeated for all photons
            vec_diff(&nodes[ea], &nodes[ec], &vecAC); //repeated for all photons
            vec_diff(&nodes[ea], &(r->p0), &vecS);
            vec_cross(&vecAB, &vecAC, &vecN); //repeated for all photons, vecN can be precomputed
            bary[facemap[i]] = -vec_dot(&vecS, &vecN);
        }

        for (i = 0; i < 4; i++) {
            if (bary[i] < -1e-4f) {
                include = 0;
            }
        }

        if (include) {
            r->eid = mesh->srcelem[is];
            float s = 0.f;

            for (i = 0; i < 4; i++) {
                s += bary[i];
            }

            r->bary0.x = bary[0] / s;
            r->bary0.y = bary[1] / s;
            r->bary0.z = bary[2] / s;
            r->bary0.w = bary[3] / s;

            for (i = 0; i < 4; i++) {
                if ((bary[i] / s) < 1e-4f) {
                    r->faceid = ifacemap[i] + 1;
                }
            }

            break;
        }
    }

    if (is == mesh->srcelemlen) {
        #pragma omp critical
        {
            MMC_FPRINTF(cfg->flog, "all tetrahedra (%d) labeled with -1 do not enclose the source!\n", mesh->srcelemlen);

            if (mesh->srcelemlen) {
                int* elems = (int*)(mesh->elem + (mesh->srcelem[0] - 1) * mesh->elemlen);
                MMC_FPRINTF(cfg->flog, "elem %d %d [%f %f %f] \n", mesh->srcelem[0], elems[0], mesh->node[elems[0] - 1].x, mesh->node[elems[0] - 1].y, mesh->node[elems[0] - 1].z);
                MMC_FPRINTF(cfg->flog, "elem %d %d [%f %f %f] \n", mesh->srcelem[0], elems[1], mesh->node[elems[1] - 1].x, mesh->node[elems[1] - 1].y, mesh->node[elems[1] - 1].z);
                MMC_FPRINTF(cfg->flog, "elem %d %d [%f %f %f] \n", mesh->srcelem[0], elems[2], mesh->node[elems[2] - 1].x, mesh->node[elems[2] - 1].y, mesh->node[elems[2] - 1].z);
                MMC_FPRINTF(cfg->flog, "elem %d %d [%f %f %f] \n", mesh->srcelem[0], elems[3], mesh->node[elems[3] - 1].x, mesh->node[elems[3] - 1].y, mesh->node[elems[3] - 1].z);
                MMC_FPRINTF(cfg->flog, "source position [%e %e %e] \n", r->p0.x, r->p0.y, r->p0.z);
                MMC_FPRINTF(cfg->flog, "bary centric volume [%e %e %e %e] \n", bary[0], bary[1], bary[2], bary[3]);
            }

            MESH_ERROR("initial element does not enclose the source!");
        }
    }
}


/**
 * @brief Deal with absorption using MCML-like algorithm (albedo-weight MC)
 *
 * This function performs MCML-like absorption calculations
 *
 * \param[in,out] r: the current ray
 * \param[in] mesh: the mesh data structure
 * \param[in,out] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

void albedoweight(ray* r, tetmesh* mesh, mcconfig* cfg, visitor* visit) {
    float* baryp0 = &(r->bary0.x);
    int eid = r->eid - 1;
    int* ee = (int*)(mesh->elem + eid * mesh->elemlen);
    int i, tshift, datalen = (cfg->method == rtBLBadouelGrid) ? cfg->crop0.z : ( (cfg->basisorder) ? mesh->nn : mesh->ne);;

    medium* prop = mesh->med + (mesh->type[eid]);
    float ww = r->weight;

    r->weight *= prop->mus / (prop->mua + prop->mus);

    ww -= r->weight;
    r->Eabsorb += ww;

    tshift = MIN( ((int)((r->photontimer - cfg->tstart) * visit->rtstep)), cfg->maxgate - 1 ) * datalen;

    if (cfg->method == rtBLBadouelGrid) {

    } else {
        if (!cfg->basisorder) {
            if (cfg->isatomic)
                #pragma omp atomic
                mesh->weight[eid + tshift] += ww;
            else {
                mesh->weight[eid + tshift] += ww;
            }
        } else {
            if (cfg->isatomic)
                for (i = 0; i < 4; i++)
                    #pragma omp atomic
                    mesh->weight[ee[i] - 1 + tshift] += ww * baryp0[i];
            else
                for (i = 0; i < 4; i++) {
                    mesh->weight[ee[i] - 1 + tshift] += ww * baryp0[i];
                }
        }
    }
}

void visitor_init(mcconfig* cfg, visitor* visit) {
    visit->launchweight = (double*)calloc(cfg->srcnum, sizeof(double));
    visit->absorbweight = (double*)calloc(cfg->srcnum, sizeof(double));
    visit->kahanc0 = (double*)calloc(cfg->srcnum, sizeof(double));
    visit->kahanc1 = (double*)calloc(cfg->srcnum, sizeof(double));
}

void visitor_clear(visitor* visit) {
    free(visit->launchweight);
    visit->launchweight = NULL;
    free(visit->absorbweight);
    visit->absorbweight = NULL;
    free(visit->kahanc0);
    visit->kahanc0 = NULL;
    free(visit->kahanc1);
    visit->kahanc1 = NULL;
}

void updateroi(int immctype, ray* r, tetmesh* mesh) {
    if (immctype == 1 && mesh->edgeroi) {
        r->roisize = (float*)(mesh->edgeroi + (r->eid - 1) * 6);
    } else if (mesh->faceroi) {
        r->roisize = (float*)(mesh->faceroi + (r->eid - 1) * 4);
    }
}