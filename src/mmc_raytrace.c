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

/***************************************************************************//**
\file    mmc_raytrace.c

\brief   Core unit for mash-based photon transport using ray tracing algorithms
*******************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mmc_raytrace.h"
#include "mmc_const.h"

/**<  Macro to enable SSE4 based ray-tracers */

#ifdef MMC_USE_SSE
    #include <smmintrin.h>
    __m128 int_coef;                            /**<  a global variable used for SSE4 ray-tracers */
#endif

#define F32N(a) ((a) & 0x80000000)          /**<  Macro to test if a floating point is negative */
#define F32P(a) ((a) ^ 0x80000000)          /**<  Macro to test if a floating point is positive */
#ifndef JUST_BELOW_ONE
    #define JUST_BELOW_ONE  0.9998f          /**<  used to clamp pattern index to avoid overflow */
#endif

static inline int mmc_is_fluence_like(char outputtype) {
    return outputtype != otEnergy && outputtype != otWP && outputtype != otWL;
}

static inline float mmc_fluence_deposit(float w0, float w1, float mua, float len) {
    return (mua < EPS) ? (w0 * len) : ((w0 - w1) / mua);
}

#ifdef MMC_USE_SSE  /* only the SSE branchless Badouel tracer makes the RF forward deposit */
/**
 * RF forward (frequency domain) deposit of one step on the CPU.  The photon's complex
 * weight is |w| e^{-i omega t}: r->weight stays the real, decaying magnitude, so the
 * roulette, escape and absorption tallies are unchanged, and the phase is omega times
 * the time of flight.  A step of length len that starts at weight w0 and phase phase0
 * deposits the fluence
 *     w0 e^{-i phase0} (1 - e^{-(mua + i a) len}) / (mua + i a),   a = omega n / c0,
 * the complex Beer-Lambert deposit the OpenCL/CUDA/OptiX kernels make.
 */
static inline void mmc_rf_deposit(float w0, double phase0, float mua, float a, float len, float* re, float* im) {
    double E = exp(-(double)mua * len), c = cos((double)a * len), s = sin((double)a * len);
    double nr = 1.0 - E * c, ni = E * s, z2 = (double)mua * mua + (double)a * a, fr, fi;

    if (z2 < (double)EPS * EPS) {
        fr = len;
        fi = 0.0;
    } else {
        fr = (nr * mua + ni * a) / z2;
        fi = (ni * mua - nr * a) / z2;
    }

    c = cos(phase0);
    s = sin(phase0);
    *re = (float)(w0 * (fr * c + fi * s));
    *im = (float)(w0 * (fi * c - fr * s));
}
#endif

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

extern const int out[4][3];     /**< defined in mmc_mesh.c, node orders for each face, in counter-clock-wise orders*/

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

inline void interppos(float3* w, float3* p1, float3* p2, float3* p3, float3* pout) {
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

inline void getinterp(float w1, float w2, float w3, FLOAT3* p1, FLOAT3* p2, FLOAT3* p3, FLOAT3* pout) {
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

void fixphoton(FLOAT3* p, FLOAT3* nodes, int* ee) {
    FLOAT3 c0 = {0.f, 0.f, 0.f};
    int i;

    /*calculate element centroid*/
    for (i = 0; i < 4; i++) {
        vec_add3(&c0, nodes + ee[i] - 1, &c0);
    }

    p->x += (c0.x * 0.25f - p->x) * FIX_PHOTON;
    p->y += (c0.y * 0.25f - p->y) * FIX_PHOTON;
    p->z += (c0.z * 0.25f - p->z) * FIX_PHOTON;
}

/**
 * \brief iMMC: trace the ROIs within the step, then take the step's medium from the state at p0
 *
 * Without an ROI crossing that is the new inroi. With one, the part before it is on the side of
 * p0, which the element's own tests (pin) decide after an element change (the stored state may
 * come from an ROI this element does not test); if the step then scatters first, the crossing
 * is dropped.
 *
 * \return the medium of the step; *mus, *rc and *dlen are updated with it
 */
static medium* roi_step(ray* r, raytracer* tracer, mcconfig* cfg, medium* prop, float Lp0, float* mus, float* rc, float* dlen) {
    int oldin = r->inroi, pin = traceroi(r, tracer, cfg->implicit, 0, cfg);

    if ((r->roitype == rtNone) ? (r->inroi != oldin) : (pin >= 0 && pin != oldin)) {
        prop = tracer->mesh->med + (((r->roitype == rtNone) ? r->inroi : pin) ? tracer->mesh->prop : tracer->mesh->type[r->eid - 1]);
        *mus = (cfg->mcmethod == mmMCX) ? prop->mus : (prop->mua + prop->mus);
        *rc = prop->n * R_C0;
        *dlen = (*mus <= EPS) ? R_MIN_MUS : r->slen / *mus;

        if (r->roitype == rtNone) {
            r->isend = (Lp0 > *dlen);
            r->Lmove = ((r->isend) ? *dlen : Lp0);
        } else if (*dlen < r->Lmove) {
            r->roitype = rtNone;
            r->inroi = pin;
            r->isend = 1;
            r->Lmove = *dlen;
        }
    }

    return prop;
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

    float3 pcrx = {0.f}, p1 = {0.f};
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
    r->lasthit = (r->roitype != rtNone && r->roiidx >= 0) ? ((r->roiidx << 3) | r->roitype) : 0;
    r->faceid = -1;
    r->isend = 0;
    r->Lmove = 0.f;
    vec_add(&(r->p0), &(r->vec), &p1);
    vec_cross(&(r->p0), &p1, &pcrx);
    ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);

    if (cfg->implicit && r->inroi) {
        prop = tracer->mesh->med + tracer->mesh->prop;
    } else {
        prop = tracer->mesh->med + (tracer->mesh->type[eid]);
    }

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

            Lp0 = dist_3((FLOAT3*) & (r->p0), &(r->pout));
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

    // implicit MMC - test if ray intersects with edge/face/node ROI boundaries
    if (cfg->implicit) {
        prop = roi_step(r, tracer, cfg, prop, Lp0, &mus, &rc, &dlen);
    }

    visit->raytet++;

    if (tracer->mesh->type[eid] == 0) {
        visit->raytet0++;
    }

    if (r->pout.x != MMC_UNDEFINED) {
        if ((int)((r->photontimer + r->Lmove * rc - cfg->tstart)*visit->rtstep) >= (int)((cfg->tend - cfg->tstart)*visit->rtstep)) { /*exit time window*/
            r->faceid = -2;
            r->pout.x = MMC_UNDEFINED;
            r->Lmove = MAX((cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f, 0.f); /* never backwards in time */
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

            if (mmc_is_fluence_like(cfg->outputtype)) {
                ww = mmc_fluence_deposit(currweight, r->weight, prop->mua, r->Lmove);
            }

            r->photontimer += r->Lmove * rc;

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * tracer->mesh->ne;
            } else {
                tshift = MIN(MAX((int)((r->photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * tracer->mesh->ne;
            }

            if (cfg->mcmethod == mmMCX && tracer->mesh->type[eid] > 0 && r->photontimer >= cfg->tstart) { /* no fluence in non-tissue (label<=0) elements, as mcx never deposits in label-0 voxels, nor before the time window */
                if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                    #pragma omp atomic
                    tracer->mesh->weight[eid + tshift] += ww;
                } else if (cfg->srctype == stPattern) { // must be pattern and srcnum more than 1
                    int pidx; // pattern index

                    for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                        #pragma omp atomic
                        tracer->mesh->weight[(eid + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
                    }
                }
            }
        } else {
            if (cfg->debuglevel & dlBary) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"barycentric\": {\"p\":[%f, %f, %f], \"eid\":%d, \"bary\":[%f, %f, %f, %f]}}\n",
                            r->photonid, r->p0.x, r->p0.y, r->p0.z, r->eid, baryout[0], baryout[1], baryout[2], baryout[3]);
            }

            if (Lp0 > EPS) {
                r->photontimer += r->Lmove * rc;
                ww = currweight - r->weight;

                if (prop->mua > 0.f || mmc_is_fluence_like(cfg->outputtype)) {
                    ratio = r->Lmove / Lp0;

                    if (prop->mua > 0.f) {
                        r->Eabsorb += ww;
                    }

                    if (mmc_is_fluence_like(cfg->outputtype)) {
                        ww = mmc_fluence_deposit(currweight, r->weight, prop->mua, r->Lmove);
                    }

                    if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                        tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * tracer->mesh->nn;
                    } else {
                        tshift = MIN(MAX((int)((r->photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * tracer->mesh->nn;
                    }

                    if (cfg->debuglevel & dlAccum) {
                        MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"add\": {\"loc\":[%f, %f, %f], \"w\":%e, \"eid\":%d, \"l\":%e}}\n",
                                    r->photonid, r->p0.x - (r->Lmove * 0.5f)*r->vec.x, r->p0.y - (r->Lmove * 0.5f)*r->vec.y, r->p0.z - (r->Lmove * 0.5f)*r->vec.z, ww, eid + 1, dlen);
                    }

                    ww *= 0.5f;

                    if (r->isend) {
                        for (i = 0; i < 4; i++) {
                            baryout[i] = (1.f - ratio) * baryp0[i] + ratio * baryout[i];
                        }
                    }

                    if (cfg->mcmethod == mmMCX && tracer->mesh->type[eid] > 0 && r->photontimer >= cfg->tstart) {
                        if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                            for (i = 0; i < 4; i++) {
                                #pragma omp atomic
                                tracer->mesh->weight[ee[i] - 1 + tshift] += ww * (baryp0[i] + baryout[i]);
                            }
                        } else if (cfg->srctype == stPattern) { // must be pattern and srcnum more than 1
                            int pidx; // pattern index

                            for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                for (i = 0; i < 4; i++) {
                                    #pragma omp atomic
                                    tracer->mesh->weight[(ee[i] - 1 + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[r->posidx * cfg->srcnum + pidx] * (baryp0[i] + baryout[i]);
                                }
                            }
                        }
                    }
                }

                if (r->isend) {
                    memcpy(baryp0, baryout, sizeof(float4));
                } else {
                    if (r->nexteid > 0 && faceidx >= 0) {
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

static inline __m128 rcp_nr(const __m128 a) {
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

static inline int havel_sse4(float3* vecN, float3* bary, const __m128 o, const __m128 d) {
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

    float3 bary = {1e10f, 0.f, 0.f, 0.f};
    MMC_ALIGN(16) float barypout[4];
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
                r->Lmove = MAX((cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f, 0.f); /* never backwards in time */
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
            }

            if (mmc_is_fluence_like(cfg->outputtype)) {
                ww = mmc_fluence_deposit(currweight, r->weight, prop->mua, r->Lmove);
            }

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            } else {
                tshift = MIN(MAX((int)((r->photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            }

            if (cfg->debuglevel & dlAccum) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"pre-add\": {\"p\":[%f, %f, %f], \"lmin\":%e, \"eid\":%d, \"l\":%e}}\n",
                            r->photonid, r->p0.x, r->p0.y, r->p0.z, bary.x, eid + 1, r->Lmove);
            }

            T = _mm_mul_ps(T, _mm_set1_ps(r->Lmove));
            S = _mm_add_ps(O, T);
            _mm_store_ps(&(r->p0.x), S);

            barypout[out[i][0]] = 1.f - bary.y - bary.z;
            barypout[out[i][1]] = bary.y;
            barypout[out[i][2]] = bary.z;
            barypout[facemap[i]] = 0.f;

            if (cfg->debuglevel & dlBary) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"barycentric\": {\"p\":[%f, %f, %f], \"eid\":%d, \"bary\":[%f, %f, %f, %f]}}\n",
                            r->photonid, r->p0.x, r->p0.y, r->p0.z, r->eid, barypout[0], barypout[1], barypout[2], barypout[3]);
            }

            T = _mm_load_ps(barypout);      /* bary centric at pout */
            O = _mm_load_ps(&(r->bary0.x)); /* bary centric at p0 */

            dlen = r->Lmove / bary.x;       /* normalized moving length */

            if (cfg->debuglevel & dlAccum) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"add\": {\"loc\":[%f, %f, %f], \"w\":%e, \"eid\":%d, \"l\":%e}}\n",
                            r->photonid, r->p0.x - (r->Lmove * 0.5f)*r->vec.x, r->p0.y - (r->Lmove * 0.5f)*r->vec.y, r->p0.z - (r->Lmove * 0.5f)*r->vec.z, ww, eid + 1, dlen);
            }

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
            if (cfg->mcmethod == mmMCX && tracer->mesh->type[eid] > 0 && r->photontimer >= cfg->tstart) {
                if (!cfg->basisorder) {
                    if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                        #pragma omp atomic
                        tracer->mesh->weight[eid + tshift] += ww;
                    } else if (cfg->srctype == stPattern) { // must be pattern and srcnum more than 1
                        int pidx; // pattern index

                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                            #pragma omp atomic
                            tracer->mesh->weight[(eid + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
                        }
                    }
                } else {
                    T = _mm_mul_ps(_mm_add_ps(O, S), _mm_set1_ps(ww * 0.5f));
                    _mm_store_ps(barypout, T);

                    if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                        for (j = 0; j < 4; j++) {
                            #pragma omp atomic
                            tracer->mesh->weight[ee[j] - 1 + tshift] += barypout[j];
                        }
                    } else if (cfg->srctype == stPattern) {
                        int pidx; // pattern index

                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                            for (j = 0; j < 4; j++) {
                                #pragma omp atomic
                                tracer->mesh->weight[(ee[j] - 1 + tshift)*cfg->srcnum + pidx] += barypout[j] * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
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

    float3 bary = {1e10f, 0.f, 0.f, 0.f};
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
            r->Lmove = MAX((cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f, 0.f); /* never backwards in time */
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

            if (mmc_is_fluence_like(cfg->outputtype)) {
                ww = mmc_fluence_deposit(currweight, r->weight, prop->mua, r->Lmove);
            }

            r->photontimer += r->Lmove * rc;

            if (cfg->outputtype == otWL || cfg->outputtype == otWP) {
                tshift = MIN( ((int)(cfg->replaytime[r->photonid] * visit->rtstep)), cfg->maxgate - 1 ) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            } else {
                tshift = MIN(MAX((int)((r->photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * (cfg->basisorder ? tracer->mesh->nn : tracer->mesh->ne);
            }

            if (cfg->debuglevel & dlAccum) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"pre-add\": {\"loc\":[%f, %f, %f], \"lmin\":%e, \"eid\":%d, \"l\":%e}}\n",
                            r->photonid, r->p0.x, r->p0.y, r->p0.z, bary.x, eid + 1, r->Lmove);
            }

            r->p0.x += r->Lmove * r->vec.x;
            r->p0.y += r->Lmove * r->vec.y;
            r->p0.z += r->Lmove * r->vec.z;

            if (cfg->mcmethod == mmMCX && tracer->mesh->type[eid] > 0 && r->photontimer >= cfg->tstart) {
                if (!cfg->basisorder) {
                    if (cfg->isatomic)
                        #pragma omp atomic
                        tracer->mesh->weight[eid + tshift] += ww;
                    else {
                        tracer->mesh->weight[eid + tshift] += ww;
                    }
                } else {
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

#endif


/*
 * Implicit MMC (iMMC) ROI tests
 *
 * Edge and node ROIs are capsules (segments with spherical caps) and spheres, listed per element
 * by mesh_roilist(); the photon is traced against their union. Face ROIs are slabs below the
 * labeled faces of an element. Each test classifies p0 as in or out by the sign of its implicit
 * function, except inside a rounding band about the surface, where the stored r->inroi (set by
 * the crossing that put the photon there) decides; the state then selects the one crossing to look
 * for (the exit if in, the entry if out). The same tests are used by the GPU (mmc_core.cl).
 */

#define ROI_QUAD_KU   (64.f * 5.9604645e-8f)   /**< rounding band of the quadric tests, per unit of |w|^2+r^2 */
#define ROI_QUAD_AMIN 1e-12f                   /**< below this, a ray is parallel to the capsule axis */
#define ROI_SLAB_KB   4.f                      /**< rounding band of the slab test, in units of roieps */
#define ROI_MAX_TINY  4                        /**< consecutive floored hits before the Fresnel event is dropped */

/**
 * \brief Solve a*t^2 + 2*b*t + c = 0 for the crossing selected by the in/out state of p0
 *
 * A positive selfeps marks the primitive crossed by the previous step, on whose surface the
 * photon sits: a root closer than selfeps is that crossing again, so the side is taken from the
 * direction of travel (b < 0: inward) and only the far root counts. An outside state (band)
 * with p0 already past the entry root and moving in is an entry missed between two steps, made at
 * t = 0.
 *
 * \return htOutIn/htInOut with *thit in [0, L], else htNoHitIn/htNoHitOut
 */
static inline int roi_quad(float a, float b, float c, float w2, float r2, float L, int inroi, float selfeps, float* thit) {
    float dc = ROI_QUAD_KU * (w2 + r2), disc = b * b - a * c, q, tn, tf, ts;
    int in = (c < -dc) ? 1 : ((c > dc) ? 0 : (inroi != 0));

    if (a < ROI_QUAD_AMIN || disc <= 0.f) {
        return (in && a < ROI_QUAD_AMIN) ? htNoHitIn : htNoHitOut;
    }

    q = -(b + copysignf(sqrtf(disc), b));
    tn = fminf(q / a, c / q);
    tf = fmaxf(q / a, c / q);
    ts = in ? tf : tn;

    if (selfeps > 0.f && ts > 0.f && ts <= L && ts < selfeps) {
        int inward = (b < 0.f);

        if (inward != (inroi != 0)) {
            *thit = 0.f;        /* the stored side is wrong: cross here */
            return inward ? htOutIn : htInOut;
        }

        if (!inward || tf < selfeps) {
            return htNoHitOut;
        }

        if (tf > L) {
            return htNoHitIn;
        }

        *thit = tf;
        return htInOut;
    }

    if (!in) {
        if (tn > 0.f && tn <= L) {
            *thit = tn;
            return htOutIn;
        }

        if (tn <= 0.f && tf > 0.f) {
            *thit = 0.f;        /* past the entry and moving in: the missed entry */
            return htOutIn;
        }

        return htNoHitOut;
    }

    if (tf <= 0.f) {
        return htNoHitOut;
    }

    if (tf <= L) {
        *thit = tf;
        return htInOut;
    }

    return htNoHitIn;
}

/** \brief The root of a*t^2 + 2*b*t + c = 0 selected by the state in, q = -(b + sign(b)*sqrt(disc)) */
#define ROI_QUAD_ROOT(in, q, a, c) ((((in) != 0) == ((q) > 0.f)) ? (q) / (a) : (c) / (q))

/** \brief The capsule root selected by the state in: the cylinder's if within the segment, else a cap's; 0 if missed */
static inline int roi_capsule_root(ray* r, FLOAT3* B, FLOAT3* w, float s0, float sv, float a, float cc, float dsc, float qc,
                                   float r2, int in, float* t) {
    int capb = ((sv > 0.f) == (in != 0));
    FLOAT3 wk;
    float bs, cs, dss;

    if (a >= ROI_QUAD_AMIN) {
        float s;

        if (dsc <= 0.f) {
            return 0;
        }

        *t = ROI_QUAD_ROOT(in, qc, a, cc);
        s = s0 + *t * sv;

        if (s >= 0.f && s <= 1.f) {
            return 1;
        }

        capb = (s > 1.f);
    }

    if (capb) {
        vec_diff3(B, (FLOAT3*)&r->p0, &wk);
    } else {
        wk = *w;
    }

    bs = wk.x * r->vec.x + wk.y * r->vec.y + wk.z * r->vec.z;
    cs = vec_dot3(&wk, &wk) - r2;
    dss = bs * bs - cs;

    if (dss <= 0.f) {
        return 0;
    }

    bs = -(bs + copysignf(sqrtf(dss), bs));
    *t = ROI_QUAD_ROOT(in, bs, 1.f, cs);
    return 1;
}

/**
 * \brief Capsule (edge ROI) test: the cylinder of segment A-B, then at most one end-cap sphere
 *
 * The cylinder root selected by the state is kept if its axial position s is in [0, 1];
 * otherwise the cap on that side (A if s < 0, B if s > 1; by the sign of v.(B-A) if the ray is
 * parallel to the axis) is solved. Self-intersections as in roi_quad().
 */
static int roi_capsule(ray* r, FLOAT3* A, FLOAT3* B, float rt, float selfeps, float* thit) {
    FLOAT3 v = {r->vec.x, r->vec.y, r->vec.z}, e, w, wp, dp, wc;
    float ie, s0, sv, sc, r2 = rt * rt, w2, dc, a, bc, cc, dsc, qc, t, L = r->Lmove;
    int in, self = 0;

    vec_diff3(A, B, &e);
    vec_diff3(A, (FLOAT3*)&r->p0, &w);
    ie = 1.f / vec_dot3(&e, &e);
    s0 = vec_dot3(&w, &e) * ie;
    sv = vec_dot3(&v, &e) * ie;
    sc = fminf(1.f, fmaxf(0.f, s0));
    vec_mult_add3(&e, &w, -s0, 1.f, &wp);       /* w and v perpendicular to the axis */
    vec_mult_add3(&e, &v, -sv, 1.f, &dp);
    vec_mult_add3(&e, &w, -sc, 1.f, &wc);       /* p0 from the closest point of the segment */
    w2 = vec_dot3(&wc, &wc);
    dc = ROI_QUAD_KU * (w2 + r2);
    in = (w2 - r2 < -dc) ? 1 : ((w2 - r2 > dc) ? 0 : (r->inroi != 0));
    a = vec_dot3(&dp, &dp);
    bc = vec_dot3(&wp, &dp);
    cc = vec_dot3(&wp, &wp) - r2;
    dsc = bc * bc - a * cc;
    qc = -(bc + copysignf(sqrtf(fmaxf(dsc, 0.f)), bc));

    for (;;) {    /* a second pass only after a self-hit from outside, for the exit root */
        int inward;

        if (!roi_capsule_root(r, B, &w, s0, sv, a, cc, dsc, qc, r2, in, &t)) {
            return htNoHitOut;
        }

        if (self || !(selfeps > 0.f && t > 0.f && t <= L && t < selfeps)) {
            break;
        }

        inward = (vec_dot3(&wc, &v) < 0.f);    /* self-intersection, see roi_quad() */

        if (inward != (r->inroi != 0)) {
            *thit = 0.f;
            return inward ? htOutIn : htInOut;
        }

        if (!inward) {
            return htNoHitOut;
        }

        self = 1;

        if (in) {
            break;
        }

        in = 1;
    }

    if (self && t <= L && t < selfeps) {
        return htNoHitOut;
    }

    if (t > 0.f && t <= L) {
        *thit = t;
        return in ? htInOut : htOutIn;
    }

    if (!in && t <= 0.f && vec_dot3(&wc, &v) < 0.f) {
        *thit = 0.f;    /* the missed entry, see roi_quad() */
        return htOutIn;
    }

    return (in && t > 0.f) ? htNoHitIn : htNoHitOut;
}

/**
 * \brief Face (slab) ROI test: the plane at depth thick below the face (point fpt, outward normal fn)
 *
 * The depth below the face along the ray is d0 - t*s, d0 = (fpt-p0).fn, s = v.fn, and the slab is
 * the depth range [0, thick]. The face plane itself is the element face, crossed by the ray-tet
 * test; a photon rounded beyond it (d0 < -band) is outside for the whole step.
 */
static inline int roi_slab(ray* r, FLOAT3* fpt, FLOAT3* fn, float thick, float band, float* thit) {
    FLOAT3 w;
    float d0, e, s, ss, es;
    int in;

    vec_diff3((FLOAT3*)&r->p0, fpt, &w);
    d0 = vec_dot3(&w, fn);

    if (d0 < -band) {
        return htNoHitOut;
    }

    e = d0 - thick;
    s = vec_dot3((FLOAT3*)&r->vec, fn);
    in = (e < -band) ? 1 : ((e > band) ? 0 : (r->inroi != 0));
    ss = in ? -s : s;
    es = in ? -e : e;

    if (ss <= 0.f || es >= r->Lmove * ss) {
        return in ? htNoHitIn : htNoHitOut;
    }

    *thit = fmaxf(e / s, 0.f);
    return in ? htInOut : htOutIn;
}

/**
 * \brief Floor an ROI hit closer than roieps so that the photon advances
 *
 * After ROI_MAX_TINY consecutive floored hits (a crease between two ROIs) the crossing is kept
 * but roiidx is cleared, so that no Fresnel event can trap the photon there.
 */
static inline void roi_tinyhit(ray* r, int hit, float roieps) {
    if (hit != htInOut && hit != htOutIn) {
        return;
    }

    if (r->Lmove < roieps) {
        r->Lmove = roieps;

        if (++r->tinyhits > ROI_MAX_TINY) {
            r->roiidx = -1;
            r->tinyhits = 0;
        }
    } else {
        r->tinyhits = 0;
    }
}

/**
 * \brief Edge/node iMMC: trace the ray against the union of the element's ROI primitives
 *
 * Each primitive is tested with its own in/out state. The photon is in the ROI at p0 if it is in
 * any primitive. An outside photon enters the union at the first entry E; an inside one leaves it
 * at the last exit X of the primitives it is in, unless it first enters another primitive
 * (E < X): that is an inner crossing (rtInner), which ends the step without a medium change.
 * r->lasthit names the primitive hit by the previous step, for the self-intersection test.
 *
 * \return the in/out state at p0 (1: in the union)
 */
static int traceroi_list(ray* r, raytracer* tracer, int doinit, mcconfig* cfg) {
    const tetmesh* mesh = tracer->mesh;
    const int beg = mesh->roioff[r->eid - 1], n = mesh->roioff[r->eid] - beg;
    const int* idx = mesh->roioff + beg;
    const int self = (!doinit && r->lasthit) ? (r->lasthit >> 3) : -1;
    float X = 0.f, E = 1e30f;
    int nin = 0, jx = -1, je = -1, k, hs;

    for (k = 0; k < n; k++) {
        float* d = mesh->roilist + (size_t)idx[k] * 8, t = 1e30f, se = (k == self) ? mesh->roieps : 0.f;
        /* p0 from the center of the bounding sphere (radius^2 in d[7]) */
        FLOAT3 w = {r->p0.x - 0.5f * (d[0] + d[4]), r->p0.y - 0.5f * (d[1] + d[5]), r->p0.z - 0.5f * (d[2] + d[6])};
        float wb = w.x * r->vec.x + w.y * r->vec.y + w.z * r->vec.z, wc = vec_dot3(&w, &w) - d[7], we = -wb - r->Lmove;
        int isin;

        if (wc > 0.f && (wb >= 0.f || wb * wb <= wc || (we > 0.f && we * we > wb * wb - wc))) {
            hs = htNoHitOut;    /* the step stays outside the bounding sphere */
        } else if (d[3] > 0.f) {
            hs = roi_capsule(r, (FLOAT3*)d, (FLOAT3*)(d + 4), d[3], se, &t);
        } else {                /* a sphere: the bounding sphere's center is its center */
            hs = roi_quad(1.f, wb, vec_dot3(&w, &w) - d[3] * d[3], vec_dot3(&w, &w), d[3] * d[3], r->Lmove, r->inroi, se, &t);
        }

        /* in this primitive: its exit (none if the step ends inside) competes for the last exit X */
        isin = (hs == htInOut || hs == htNoHitIn);
        t = (hs == htNoHitIn) ? 1e30f : t;
        nin += isin;

        if (isin && (jx < 0 || t > X)) {
            X = t, jx = k;
        }

        if (hs == htOutIn && t < E) {
            E = t, je = k;
        }
    }

    r->roitype = rtNone;

    if (doinit) {
        r->inroi |= (nin > 0);
        return (nin > 0);
    }

    hs = htNone;

    if (je >= 0 && (nin == 0 || E < X)) {
        hs = htOutIn, k = je, r->Lmove = E;
    } else if (nin > 0 && X < 1e30f) {
        hs = htInOut, k = jx, r->Lmove = X;
    }

    if (hs != htNone) {
        r->inroi = (hs == htOutIn);
        r->roiidx = k;
        r->roitype = (nin > 0 && hs == htOutIn) ? rtInner : ((mesh->roilist[(size_t)idx[k] * 8 + 3] > 0.f) ? rtEdge : rtNode);
    } else {
        r->inroi = (nin > 0);
    }

    roi_tinyhit(r, hs, mesh->roieps);
    return (nin > 0);
}

/**
 * \brief Implicit MMC ray-tracing core function
 *
 * With doinit = 1, only initializes r->inroi at the launch position; otherwise shortens
 * r->Lmove to the first ROI crossing, if any, and updates r->{inroi, roitype, roiidx}.
 *
 * \return edge/node mode: the in/out state at p0 from the element's own tests, which may differ
 *         from the stored one after entering an element; -1 in face mode
 */
int traceroi(ray* r, raytracer* tracer, int roitype, int doinit, mcconfig* cfg) {
    int eid = r->eid - 1, pin = -1;

    if (roitype == 1) {
        pin = traceroi_list(r, tracer, doinit, cfg);
    } else if (tracer->mesh->faceroi) {
        int* ee = (int*)(tracer->mesh->elem + eid * tracer->mesh->elemlen);
        int hs, firsthit = htNone, anyin = 0, nohitin = 0, tested = 0, i;
        float minL = r->Lmove, band = ROI_SLAB_KB * tracer->mesh->roieps;

        for (i = 0; i < 4; i++) {
            if (r->roisize[i] > 0.f) {
                float* nv = &(tracer->n[eid << 2].x), lr = r->Lmove;
                FLOAT3 fn = {nv[ifaceorder[i]], nv[ifaceorder[i] + 4], nv[ifaceorder[i] + 8]};

                hs = roi_slab(r, tracer->mesh->node + ee[nc[ifaceorder[i]][0]] - 1, &fn, r->roisize[i], band, &lr);
                tested = 1;
                anyin |= (hs == htInOut || hs == htNoHitIn);
                nohitin |= (hs == htNoHitIn);

                /* the nearest crossing wins; without one the photon is in if it is in any slab */
                if (!doinit && (hs == htInOut || hs == htOutIn) && lr < minL) {
                    minL = lr, firsthit = hs, r->roiidx = i;
                }
            }
        }

        if (doinit) {
            r->inroi |= anyin;
        } else {
            r->Lmove = minL;
            r->inroi = (firsthit != htNone) ? (firsthit == htOutIn) : (tested ? nohitin : r->inroi);
        }

        r->roitype = (firsthit != htNone) ? rtFace : rtNone;
        roi_tinyhit(r, firsthit, tracer->mesh->roieps);
    }

    if (cfg->debuglevel & dlMove) {
        MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"implicit\": {\"p\":[%f, %f, %f], \"v\":[%f, %f, %f], \"l\":%f, \"inroi\":%d, \"roitype\":%d, \"roiidx\":%d}}\n",
                    r->photonid, r->p0.x, r->p0.y, r->p0.z, r->vec.x, r->vec.y, r->vec.z, r->Lmove, r->inroi, r->roitype, r->roiidx);
    }

    return pin;
}

/**
 * \brief Face iMMC: classify the crossing of the mesh face oldeid -> eid (faceid of oldeid)
 *
 * A face labeled from one side only is the outer surface of a single-sided slab, an interface
 * between the background and the ROI medium; a face labeled from both sides is inside the slab.
 *
 * \return 1 if a photon outside the ROI enters the slab through its labeled face, 2 if a photon in
 *         the ROI leaves it there, 0 otherwise
 */
static int slabfacecross(tetmesh* mesh, int oldeid, int eid, int faceid, int inroi) {
    int* nb = mesh->facenb + (eid - 1) * mesh->elemlen;
    int lold = (mesh->faceroi[((oldeid - 1) << 2) + faceid] > 0.f), lnew = 0, k;

    for (k = 0; k < 4; k++) {
        lnew |= (nb[k] == oldeid) && (mesh->faceroi[((eid - 1) << 2) + k] > 0.f);
    }

    return (!inroi && lnew && !lold) ? 1 : ((inroi && lold && !lnew) ? 2 : 0);
}

#ifdef MMC_USE_SSE

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

    float3 bary = {1e10f, 0.f, 0.f, 0.f};
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
    r->lasthit = (r->roitype != rtNone && r->roiidx >= 0) ? ((r->roiidx << 3) | r->roitype) : 0;
    r->roitype = rtNone;
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

    /* face iMMC: reset inroi from the depths below the faces (faceroi column k: T[ifaceorder[k]]):
       in if decisively inside a slab, out if decisively outside all, else the stored state */
    if (cfg->implicit == 2 && r->roisize && r->roisize[0] != 0.f) {    /* column 0 is nonzero in an element with slabs */
        float MMC_ALIGN(16) Td[4], bd = ROI_SLAB_KB * tracer->mesh->roieps, *th = r->roisize;
        int rin = 0, rout = 1, k;

        _mm_store_ps(Td, T);

        for (k = 0; k < 4; k++) {
            float tk = Td[ifaceorder[k]];
            rin |= (tk >= -bd) & (tk < th[k] - bd);
            rout &= (th[k] <= 0.f) | (tk < -bd) | (tk > th[k] + bd);
        }

        r->inroi = rin | (r->inroi & !rout);
    }

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
        int pidx; // pattern index

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
            prop = roi_step(r, tracer, cfg, prop, Lp0, &mus, &rc, &dlen);
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
            r->Lmove = MAX((cfg->tend - r->photontimer) / (prop->n * R_C0) - 1e-4f, 0.f); /* never backwards in time */
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
                tshift = MIN(MAX((int)((r->photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * framelen;
            }

            if (cfg->debuglevel & dlAccum) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"pre-add\": {\"p\":[%f, %f, %f], \"lmin\":%e, \"eid\":%d, \"l\":%e}}\n",
                            r->photonid, r->p0.x, r->p0.y, r->p0.z, bary.x, eid + 1, r->Lmove);
            }

            if (prop->mua > 0.f) {
                r->Eabsorb += ww;
            }

            if (mmc_is_fluence_like(cfg->outputtype)) {
                ww = mmc_fluence_deposit(currweight, r->weight, prop->mua, r->Lmove);
            }

            /* RF forward: see mmc_rf_deposit; the step started at photontimer - Lmove*rc */
            int isrf = (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE && cfg->exportadjoint != NULL);
            float ww_im = 0.f, rfa = cfg->omega * rc;
            double rfphase0 = (double)cfg->omega * (r->photontimer - r->Lmove * rc);

            if (isrf && cfg->method == rtBLBadouel) {
                mmc_rf_deposit(currweight, rfphase0, prop->mua, rfa, r->Lmove, &ww, &ww_im);
            }

            T = _mm_set1_ps(r->Lmove);
            T = _mm_add_ps(S, _mm_mul_ps(O, T));
            _mm_store_ps(&(r->p0.x), T);

            if (cfg->mcmethod == mmMCX && tracer->mesh->type[eid] > 0 && r->photontimer >= cfg->tstart) {
                if (!cfg->basisorder) {
                    if (cfg->method == rtBLBadouel) {
                        unsigned int newidx = eid + tshift;
                        r->oldidx = (r->oldidx == 0xFFFFFFFF) ? newidx : r->oldidx;

                        if (newidx != r->oldidx) {
                            if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                                #pragma omp atomic
                                tracer->mesh->weight[r->oldidx] += r->oldweight;

                                if (isrf) {
                                    #pragma omp atomic
                                    cfg->exportadjoint[r->oldidx] += (float)r->oldweight_im;
                                }
                            } else if (cfg->srctype == stPattern) {
                                for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                    #pragma omp atomic
                                    tracer->mesh->weight[r->oldidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
                                }
                            }

                            r->oldidx = newidx;
                            r->oldweight = ww;
                            r->oldweight_im = ww_im;
                        } else {
                            r->oldweight += ww;
                            r->oldweight_im += ww_im;
                        }

                        if (r->faceid == -2 || !r->isend) {
                            if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                                #pragma omp atomic
                                tracer->mesh->weight[newidx] += r->oldweight;

                                if (isrf) {
                                    #pragma omp atomic
                                    cfg->exportadjoint[newidx] += (float)r->oldweight_im;
                                }
                            } else if (cfg->srctype == stPattern) {
                                for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                    #pragma omp atomic
                                    tracer->mesh->weight[newidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
                                }
                            }

                            r->oldweight = 0.f;
                            r->oldweight_im = 0.f;
                        }
                    } else {
                        float dstep, segloss, w0;
                        MMC_ALIGN(16) int4 idx;
                        int i, seg = (int)(r->Lmove / cfg->steps.x) + 1;
                        seg = (seg << 1);
                        dstep = r->Lmove / seg;
                        segloss = expf(-prop->mua * dstep);
                        /* RF: each sub-segment deposits from its own starting weight and phase */
                        float rfds = dstep, rfw = currweight, segre = 0.f, segim = 0.f;
                        double rfph = rfphase0;
                        T =  _mm_mul_ps(O, _mm_set1_ps(dstep)); /*step*/
                        O =  _mm_sub_ps(S, _mm_load_ps(&(tracer->mesh->nmin.x)));
                        S =  _mm_add_ps(O, _mm_mul_ps(T, _mm_set1_ps(0.5f))); /*starting point*/
                        dstep = 1.f / cfg->steps.x;
                        totalloss = (totalloss == 0.f) ? 0.f : (1.f - segloss) / totalloss;
                        w0 = ww;

                        for (i = 0; i < seg; i++) {
                            P = _mm_cvttps_epi32(_mm_mul_ps(S, _mm_set1_ps(dstep)));
                            _mm_store_si128((__m128i*) & (idx.x), P);

                            /* a midpoint rounded just outside the grid, as an extreme refractive
                               index or the rounding of the steps can make it, must not index past
                               it; one unsigned test per axis, the clamp itself is rare */
                            if ((unsigned)idx.x >= cfg->dim.x || (unsigned)idx.y >= cfg->dim.y || (unsigned)idx.z >= cfg->dim.z) {
                                idx.x = MIN(MAX(idx.x, 0), (int)cfg->dim.x - 1);
                                idx.y = MIN(MAX(idx.y, 0), (int)cfg->dim.y - 1);
                                idx.z = MIN(MAX(idx.z, 0), (int)cfg->dim.z - 1);
                            }

                            unsigned int newidx = idx.z * cfg->crop0.y + idx.y * cfg->crop0.x + idx.x + tshift;
                            r->oldidx = (r->oldidx == 0xFFFFFFFF) ? newidx : r->oldidx;

                            if (isrf) {
                                mmc_rf_deposit(rfw, rfph, prop->mua, rfa, rfds, &segre, &segim);
                                rfw *= segloss;
                                rfph += (double)rfa * rfds;
                            }

                            if (newidx != r->oldidx) {
                                if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                                    #pragma omp atomic
                                    tracer->mesh->weight[r->oldidx] += r->oldweight;

                                    if (isrf) {
                                        #pragma omp atomic
                                        cfg->exportadjoint[r->oldidx] += (float)r->oldweight_im;
                                    }
                                } else if (cfg->srctype == stPattern) {
                                    for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                            #pragma omp atomic
                                            tracer->mesh->weight[r->oldidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
                                        }
                                    }
                                }

                                r->oldidx = newidx;
                                r->oldweight = isrf ? segre : w0 * totalloss;
                                r->oldweight_im = segim;
                            } else {
                                r->oldweight += isrf ? segre : w0 * totalloss;
                                r->oldweight_im += segim;
                            }

                            if (r->faceid == -2 || !r->isend) {
                                if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                                    #pragma omp atomic
                                    tracer->mesh->weight[newidx] += r->oldweight;

                                    if (isrf) {
                                        #pragma omp atomic
                                        cfg->exportadjoint[newidx] += (float)r->oldweight_im;
                                    }
                                } else if (cfg->srctype == stPattern) {
                                    for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                                            #pragma omp atomic
                                            tracer->mesh->weight[newidx * cfg->srcnum + pidx] += r->oldweight * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
                                        }
                                    }
                                }

                                r->oldweight = 0.f;
                                r->oldweight_im = 0.f;
                            }

                            w0 *= segloss;
                            S = _mm_add_ps(S, T);
                        }
                    }
                } else {
                    int i;
                    ww *= 1.f / 3.f;
                    ww_im *= 1.f / 3.f;

                    if (cfg->srctype != stPattern || cfg->srcnum == 1) {
                        for (i = 0; i < 3; i++) {
                            #pragma omp atomic
                            tracer->mesh->weight[ee[out[faceidx][i]] - 1 + tshift] += ww;

                            if (isrf) {
                                #pragma omp atomic
                                cfg->exportadjoint[ee[out[faceidx][i]] - 1 + tshift] += ww_im;
                            }
                        }
                    } else if (cfg->srctype == stPattern) {
                        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
                            for (i = 0; i < 3; i++) {
                                #pragma omp atomic
                                tracer->mesh->weight[(ee[out[faceidx][i]] - 1 + tshift)*cfg->srcnum + pidx] += ww * cfg->srcpattern[r->posidx * cfg->srcnum + pidx];
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
 * @brief Saving photon trajectory data for debugging purposes
 * @param[in] p: the position/weight of the current photon packet
 * @param[in] id: the global index of the photon
 * @param[in] gdebugdata: pointer to the global-memory buffer to store the trajectory info
 */

void savedebugdata(ray* r, unsigned int id, mcconfig* cfg) {
    unsigned int pos;
    float* gdebugdata = cfg->exportdebugdata;

#ifdef _MSC_VER                 // MSVC knows no "atomic capture" (OpenMP 3.1)
    #pragma omp critical
#else
    #pragma omp atomic capture
#endif
    pos = cfg->debugdatalen++;

    if (pos < cfg->maxjumpdebug) {
        pos *= MCX_DEBUG_REC_LEN;
        ((unsigned int*)gdebugdata)[pos++] = id;
        gdebugdata[pos++] = r->p0.x;
        gdebugdata[pos++] = r->p0.y;
        gdebugdata[pos++] = r->p0.z;
        gdebugdata[pos++] = r->weight;
        ((unsigned int*)gdebugdata)[pos++] = r->eid;
    }
}

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

    if (cfg->debuglevel & dlTraj) {
        savedebugdata(&r, (unsigned int)id, cfg);
    }

    /*use Kahan summation to accumulate weight, otherwise, counter stops at 16777216*/
    /*http://stackoverflow.com/questions/2148149/how-to-sum-a-large-number-of-float-number*/
    int pidx;

    if (cfg->srctype != stPattern || cfg->srcnum == 1) {
        r.partialpath[visit->reclen - 2] = r.weight;

        if (cfg->seed == SEED_FROM_FILE && (cfg->outputtype == otWL || cfg->outputtype == otWP)) {
            kahany = cfg->replayweight[r.photonid] - visit->kahanc0[0];    /* when replay mode, accumulate detected photon weight */
        } else {
            kahany = r.weight - visit->kahanc0[0];
        }

        kahant = visit->launchweight[0] + kahany;
        visit->kahanc0[0] = (kahant - visit->launchweight[0]) - kahany;
        visit->launchweight[0] = kahant;
    } else if (cfg->srctype == stPattern) {
        *((int*)(r.partialpath + visit->reclen - 2)) = r.posidx;

        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
            if (cfg->seed == SEED_FROM_FILE && (cfg->outputtype == otWL || cfg->outputtype == otWP)) {
                kahany = cfg->replayweight[r.photonid] - visit->kahanc0[pidx];    /* when replay mode, accumulate detected photon weight */
            } else {
                kahany = r.weight * cfg->srcpattern[r.posidx * cfg->srcnum + pidx] - visit->kahanc0[pidx];
            }

            kahant = visit->launchweight[pidx] + kahany;
            visit->kahanc0[pidx] = (kahant - visit->launchweight[pidx]) - kahany;
            visit->launchweight[pidx] = kahant;
        }
    }

#ifdef MMC_USE_SSE
    const float int_coef_arr[4] = { -1.f, -1.f, -1.f, 1.f };
    int_coef = _mm_load_ps(int_coef_arr);
#endif

    /** retrieve the iMMC ROI size and ray location at initial launch */
    if (cfg->implicit) {
        updateroi(cfg->implicit, &r, tracer->mesh);
        traceroi(&r, tracer, cfg->implicit, 1, cfg);
    }

    while (1) { /*propagate a photon until exit*/
        if (cfg->implicit) {
            updateroi(cfg->implicit, &r, tracer->mesh);
        }

        r.slen = (*tracercore)(&r, tracer, cfg, visit);

        if (r.pout.x == MMC_UNDEFINED) {
            if (r.faceid == -2) {
                break;    /*reaches the time limit; the tracer leaves pout undefined*/
            }

            /* fixphoton nudges the photon towards the centroid of its element;
               one the tracer has moved out of the mesh (eid <= 0) or lost has none */
            if (r.eid > 0 && r.eid <= (int)mesh->ne && fixcount++ < MAX_TRIAL) {
                fixphoton((FLOAT3*)&r.p0, mesh->node, (int*)(mesh->elem + (r.eid - 1)*mesh->elemlen));
                continue;
            }

            r.eid = ID_UNDEFINED;
            r.faceid = -1;
        }

        if (cfg->issavedet && r.Lmove > 0.f && r.faceid >= 0 && mesh->type[r.eid - 1] > 0) { /* faceid first: a photon the tracer gave up on has faceid -1 and eid ID_UNDEFINED */
            r.partialpath[mesh->prop - 1 + ((cfg->implicit && ROI_STEP_IN(r.inroi, r.roitype)) ? mesh->prop : mesh->type[r.eid - 1])] += r.Lmove;    /*second medianum block is the partial path*/

            if (cfg->debuglevel & dlMove) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"partialpath-1\": {\"l\":%f, \"inroi\":%d, \"roitype\":%d}}\n",
                            r.photonid, r.Lmove, r.inroi, r.roitype);
            }
        }

        if (r.faceid == -2) {
            break;    /*reaches the time limit*/
        }

        if (cfg->implicit && cfg->isreflect && r.roitype && r.roitype != rtInner && r.roiidx >= 0 && r.eid != ID_UNDEFINED && (mesh->med[cfg->his.maxmedia].n != mesh->med[mesh->type[r.eid - 1]].n)) {
            if (cfg->debuglevel & dlTraj) {
                savedebugdata(&r, (unsigned int)id, cfg);
            }

            reflectrayroi(cfg, (FLOAT3*)&r.vec, (FLOAT3*)&r.p0, tracer, &r.eid, &r.inroi, ran, r.roitype, r.roiidx);
            vec_mult_add(&r.p0, &r.vec, 1.0f, mesh->roieps, &r.p0);
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
                /* face iMMC: the labeled face of a single-sided slab separates the background from the
                   ROI medium (sf = 1: entering the slab, 2: leaving it) */
                int sf = (cfg->implicit == 2 && mesh->nslab1 && r.eid > 0) ? slabfacecross(mesh, oldeid, r.eid, r.faceid, r.inroi) : 0;

                if (sf) {
                    float nbg = mesh->med[mesh->type[(sf == 1 ? oldeid : r.eid) - 1]].n, nroi = mesh->med[mesh->prop].n;

                    if (cfg->isreflect && nbg != nroi) {
                        reflectray_n(cfg, &r.vec, tracer, &oldeid, &r.eid, r.faceid, ran, r.inroi, (sf == 1) ? nbg : nroi, (sf == 1) ? nroi : nbg);
                    }

                    if (r.eid != oldeid) {    /* transmitted */
                        r.inroi = (sf == 1);
                    }
                } else if (cfg->isreflect && (r.eid <= 0 || (!r.inroi && mesh->med[mesh->type[r.eid - 1]].n != mesh->med[mesh->type[oldeid - 1]].n))) {
                    /* in the ROI, the photon is in the ROI medium on both sides of an interior face */
                    float nold = r.inroi ? mesh->med[mesh->prop].n : mesh->med[mesh->type[oldeid - 1]].n;

                    if (! (r.eid <= 0 && ((nold == cfg->nout && cfg->isreflect != (int)bcMirror) || cfg->isreflect == (int)bcAbsorbExterior) ) ) {
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
                    MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"enter\": {\"p\":[%f, %f, %f], \"v\":[%f, %f, %f], \"w\":%f, \"eid\":%d}}\n", id, r.p0.x, r.p0.y, r.p0.z,
                                r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid);

                if (!cfg->voidtime) {
                    r.photontimer = 0.f;
                }
            }

            /*when a photon exits the domain into the background*/
            if (mesh->type[oldeid - 1] && mesh->type[r.eid - 1] == 0) {
                if (cfg->debuglevel & dlExit)
                    MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"leave\": {\"p\":[%f, %f, %f], \"v\":[%f, %f, %f], \"w\":%f, \"eid\":%d}}\n", id, r.p0.x, r.p0.y, r.p0.z,
                                r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid);

                if (!cfg->isextdet) {
                    r.eid = 0;
                    break;
                }
            }

            /*print intersection position*/
            if (r.pout.x != MMC_UNDEFINED && (cfg->debuglevel & dlMove)) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"hit\": {\"pout\":[%f, %f, %f], \"eid\":%d, \"slen\":%e}}\n", id, r.pout.x, r.pout.y, r.pout.z, r.eid, r.slen);
            }

            if (cfg->implicit) {
                updateroi(cfg->implicit, &r, tracer->mesh);
            }

            r.slen = (*tracercore)(&r, tracer, cfg, visit);

            if (cfg->issavedet && r.Lmove > 0.f && mesh->type[r.eid - 1] > 0) {
                r.partialpath[mesh->prop - 1 + ((cfg->implicit && ROI_STEP_IN(r.inroi, r.roitype)) ? mesh->prop : mesh->type[r.eid - 1])] += r.Lmove;

                if (cfg->debuglevel & dlMove) {
                    MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"partialpath-2\": {\"l\":%f, \"inroi\":%d, \"roitype\":%d}}\n",
                                r.photonid, r.Lmove, r.inroi, r.roitype);
                }
            }

            fixcount = 0;

            while (r.pout.x == MMC_UNDEFINED && r.eid > 0 && r.eid <= (int)mesh->ne && fixcount++ < MAX_TRIAL) {
                fixphoton((FLOAT3*)&r.p0, mesh->node, (int*)(mesh->elem + (r.eid - 1)*mesh->elemlen));
                r.slen = (*tracercore)(&r, tracer, cfg, visit);

                if (cfg->issavedet && r.Lmove > 0.f && mesh->type[r.eid - 1] > 0) {
                    r.partialpath[mesh->prop - 1 + ((cfg->implicit && ROI_STEP_IN(r.inroi, r.roitype)) ? mesh->prop : mesh->type[r.eid - 1])] += r.Lmove;

                    if (cfg->debuglevel & dlMove) {
                        MMC_FPRINTF(cfg->flog, "{\"id\":%u, \"partialpath-fixphoton\": {\"l\":%f, \"inroi\":%d, \"roitype\":%d}}\n",
                                    r.photonid, r.Lmove, r.inroi, r.roitype);
                    }
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
                MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"exit-boundary\": {\"p\":[%f, %f, %f], \"eid\":%d, \"slen\":%e}}\n", id, r.p0.x, r.p0.y, r.p0.z, r.eid, r.slen);
            }

            if (r.eid != ID_UNDEFINED) {
                if (cfg->debuglevel & dlExit)
                    MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"terminate\": {\"p\":[%f, %f, %f], \"v\":[%f, %f, %f], \"w\":%f, \"eid\":%d}}\n", id, r.p0.x, r.p0.y, r.p0.z,
                                r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid);

                if (cfg->issavedet && cfg->issaveexit) {                                   /*when issaveexit is set to 1*/
                    memcpy(r.partialpath + (visit->reclen - 2 - 6), &(r.p0.x), sizeof(float) * 3); /*columns 7-5 from the right store the exit positions*/
                    memcpy(r.partialpath + (visit->reclen - 2 - 3), &(r.vec.x), sizeof(float) * 3); /*columns 4-2 from the right store the exit dirs*/
                }

                if (cfg->issaveref && r.eid < 0 && mesh->dref) {
                    int tshift = MIN(MAX((int)((r.photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * mesh->nf;
                    #pragma omp atomic
                    mesh->dref[((-r.eid) - 1) + tshift] += r.weight;
                }
            } else if (r.faceid == -2 && (cfg->debuglevel & dlMove)) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"terminate-by-time\": {\"p\":[%f, %f, %f], \"eid\":%d, \"slen\":%e}}\n", id, r.p0.x, r.p0.y, r.p0.z, r.eid, r.slen);
            } else if (r.eid && r.faceid != -2  && cfg->debuglevel & dlEdge) {
                MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"terminate-by-exit\": {\"p\":[%f, %f, %f], \"eid\":%d, \"slen\":%e}}\n", id, r.p0.x, r.p0.y, r.p0.z, r.eid, r.slen);
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
            MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"move\": {\"p\":[%f, %f, %f], \"eid\":%d, \"slen\":%e, \"l\":%f}}\n", id, r.p0.x, r.p0.y, r.p0.z, r.eid, r.slen, r.Lmove);
        }

        if (cfg->minenergy > 0.f && r.weight < cfg->minenergy && (cfg->tend - cfg->tstart)*visit->rtstep <= 1.f) { /*Russian Roulette*/
            if (rand_do_roulette(ran)*cfg->roulettesize <= 1.f) {
                r.weight *= cfg->roulettesize;

                if (cfg->debuglevel & dlWeight) {
                    MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"roulette-live\": {\"w\":%e}}\n", id, r.weight);
                }
            } else {
                if (cfg->debuglevel & dlWeight) {
                    MMC_FPRINTF(cfg->flog, "{\"id\":%zu, \"roulette-die\": {\"w\":%e}}\n", id, r.weight);
                }

                break;
            }
        }

        if (cfg->implicit && cfg->isreflect && r.roitype && r.roitype != rtInner && r.roiidx >= 0 && r.eid != ID_UNDEFINED && (mesh->med[cfg->his.maxmedia].n != mesh->med[mesh->type[r.eid - 1]].n)) {
            if (cfg->debuglevel & dlTraj) {
                savedebugdata(&r, (unsigned int)id, cfg);
            }

            reflectrayroi(cfg, (FLOAT3*)&r.vec, (FLOAT3*)&r.p0, tracer, &r.eid, &r.inroi, ran, r.roitype, r.roiidx);
            vec_mult_add(&r.p0, &r.vec, 1.0f, mesh->roieps, &r.p0);
            continue;
        } else if (cfg->implicit && r.roitype) {
            continue;
        }

        mom = 0.f;
        r.slen0 = mc_next_scatter(mesh->med[(cfg->implicit && r.inroi) ? mesh->prop : mesh->type[r.eid - 1]].g, &r.vec, ran, ran0, cfg, &mom);
        r.slen = r.slen0;
        r.tinyhits = 0;

        if (cfg->debuglevel & dlTraj) {
            savedebugdata(&r, (unsigned int)id, cfg);
        }

        if (cfg->mcmethod != mmMCX) {
            albedoweight(&r, mesh, cfg, visit);
        }

        if (cfg->ismomentum && mesh->type[r.eid - 1] > 0) {              /*when ismomentum is set to 1*/
            r.partialpath[(mesh->prop << 1) - 1 + (r.inroi ? mesh->prop : mesh->type[r.eid - 1])] += mom;    /*the third medianum block stores the momentum transfer*/
        }

        /* a tissue element only, as for the momentum above: a label of 0 or
           less has no slot, and a background given a scattering coefficient
           scatters photons there, which indexed one slot before the array */
        if (r.inroi || mesh->type[r.eid - 1] > 0) {
            r.partialpath[(r.inroi ? mesh->prop : mesh->type[r.eid - 1]) - 1]++;                      /*the first medianum block stores the scattering event counts*/
        }
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
            memcpy((char*)visit->photonseed + visit->bufpos * (sizeof(RandType)*RAND_BUF_LEN), r.photonseed, (sizeof(RandType)*RAND_BUF_LEN));
        }

        visit->bufpos++;
    }

    free(r.partialpath);

    if (r.photonseed) {
        free(r.photonseed);
    }

    if (cfg->debuglevel & dlTraj) {
        savedebugdata(&r, (unsigned int)id, cfg);
    }

    if (cfg->srctype != stPattern || cfg->srcnum == 1) {
        kahany = r.Eabsorb - visit->kahanc1[0];
        kahant = visit->absorbweight[0] + kahany;
        visit->kahanc1[0] = (kahant - visit->absorbweight[0]) - kahany;
        visit->absorbweight[0] = kahant;
    } else if (cfg->srctype == stPattern) {
        for (pidx = 0; pidx < cfg->srcnum; pidx++) {
            kahany = r.Eabsorb * cfg->srcpattern[r.posidx * cfg->srcnum + pidx] - visit->kahanc1[pidx];
            kahant = visit->absorbweight[pidx] + kahany;
            visit->kahanc1[pidx] = (kahant - visit->absorbweight[pidx]) - kahany;
            visit->absorbweight[pidx] = kahant;
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

float reflectrayroi(mcconfig* cfg, FLOAT3* c0, FLOAT3* ph, raytracer* tracer, int* eid, int* inroi, RandType* ran, int roitype, int roiidx) {

    /*to handle refractive index mismatch*/
    FLOAT3 pnorm = {0.f}, *pn = &pnorm;
    float Icos, Re, Im, Rtotal, tmp0, tmp1, tmp2, n1, n2;

    if (roitype == rtEdge || roitype == rtNode) { // the outward normal of the capsule or sphere at ph
        const int* ro = tracer->mesh->roioff;
        float* d = tracer->mesh->roilist + (size_t)ro[ro[*eid - 1] + roiidx] * 8;
        FLOAT3 u;

        vec_diff3((FLOAT3*)d, ph, pn);

        if (roitype == rtEdge) {    // from the closest point of the segment (side or end cap)
            vec_diff3((FLOAT3*)d, (FLOAT3*)(d + 4), &u);
            tmp0 = fminf(1.f, fmaxf(0.f, vec_dot3(pn, &u) / vec_dot3(&u, &u)));
            vec_mult_add3(&u, pn, -tmp0, 1.f, pn);
        }

        tmp0 = 1.f / sqrtf(vec_dot3(pn, pn));
        vec_mult3(pn, tmp0, pn);
    } else {         // hit face roi: faceroi column roiidx is the face ifaceorder[roiidx]
        int baseid = (*eid - 1) << 2;

        pn->x = (&(tracer->n[baseid].x))[ifaceorder[roiidx]];
        pn->y = (&(tracer->n[baseid].x))[ifaceorder[roiidx] + 4];
        pn->z = (&(tracer->n[baseid].x))[ifaceorder[roiidx] + 8];
    }

    /*pn pointing outward*/

    // /*compute the cos of the incidence angle*/
    Icos = fabs(vec_dot3(c0, pn));

    if (*inroi) { // out->in
        n1 = tracer->mesh->med[tracer->mesh->type[*eid - 1]].n;
        n2 = tracer->mesh->med[cfg->his.maxmedia].n;

        if (roitype == rtEdge || roitype == rtNode) {
            vec_mult3(pn, -1.f, pn);
        }
    } else {     // in->out
        n1 = tracer->mesh->med[cfg->his.maxmedia].n;
        n2 = tracer->mesh->med[tracer->mesh->type[*eid - 1]].n;

        if (roitype == rtFace) {
            vec_mult3(pn, -1.f, pn);
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
            vec_mult_add3(pn, c0, -2.f * Icos, 1.f, c0);
            *inroi = !(*inroi);

            if (cfg->debuglevel & dlReflect) {
                MMC_FPRINTF(cfg->flog, "{\"implicit-reflect\": {\"norm\":[%f, %f, %f], \"v\":[%f, %f, %f], \"eid\":%d, \"inroi\":%d, \"roitype\":%d, \"roiidx\":%d}}\n", pn->x, pn->y, pn->z, c0->x, c0->y, c0->z, *eid, *inroi, roitype, roiidx);
            }
        } else if (cfg->isspecular == 2 && *eid == 0) {
            // if do transmission, but next neighbor is 0, terminate
        } else {                             /*do transmission*/
            vec_mult_add3(pn, c0, -Icos, 1.f, c0);
            vec_mult_add3(pn, c0, tmp2, n1 / n2, c0);

            if (cfg->debuglevel & dlReflect) {
                MMC_FPRINTF(cfg->flog, "{\"implicit-transmit\": {\"norm\":[%f, %f, %f], \"v\":[%f, %f, %f], \"eid\":%d, \"inroi\":%d, \"roitype\":%d, \"roiidx\":%d}}\n", pn->x, pn->y, pn->z, c0->x, c0->y, c0->z, *eid, *inroi, roitype, roiidx);
            }
        }
    } else { /*total internal reflection*/
        vec_mult_add3(pn, c0, -2.f * Icos, 1.f, c0);
        *inroi = !(*inroi);

        if (cfg->debuglevel & dlReflect) {
            MMC_FPRINTF(cfg->flog, "{\"implicit-totalreflect\": {\"norm\":[%f, %f, %f], \"v\":[%f, %f, %f], \"eid\":%d, \"inroi\":%d, \"roitype\":%d, \"roiidx\":%d}}\n", pn->x, pn->y, pn->z, c0->x, c0->y, c0->z, *eid, *inroi, roitype, roiidx);
        }
    }

    tmp0 = 1.f / sqrtf(vec_dot3(c0, c0));
    vec_mult3(c0, tmp0, c0);
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

float reflectray(mcconfig* cfg, float3* c0, raytracer* tracer, int* oldeid, int* eid, int faceid, RandType* ran, int inroi) {
    float n1, n2;

    if (cfg->implicit && inroi != 0) {    /* in the ROI medium; an ROI touching the domain surface faces nout */
        n1 = tracer->mesh->med[cfg->his.maxmedia].n;
        n2 = (*eid > 0) ? n1 : cfg->nout;
    } else {
        n1 = (*oldeid != *eid) ? tracer->mesh->med[tracer->mesh->type[*oldeid - 1]].n : cfg->nout;
        n2 = (*eid > 0) ? tracer->mesh->med[tracer->mesh->type[*eid - 1]].n : cfg->nout;
    }

    return reflectray_n(cfg, c0, tracer, oldeid, eid, faceid, ran, inroi, n1, n2);
}

/**
 * @brief Reflection/transmission at a mesh face between refractive indices n1 (oldeid side) and n2
 *
 * reflectray() with the indices given, e.g. for the labeled face of a face-iMMC slab
 */

float reflectray_n(mcconfig* cfg, float3* c0, raytracer* tracer, int* oldeid, int* eid, int faceid, RandType* ran, int inroi, float n1, float n2) {
    /*to handle refractive index mismatch*/
    float3 pnorm = {0.f}, *pn = &pnorm;
    float Icos, Re, Im, Rtotal, tmp0, tmp1, tmp2;
    int offs = (*oldeid - 1) << 2;

    faceid = ifaceorder[faceid];

    /*calculate the normal direction of the intersecting triangle*/
    if (cfg->method == rtPlucker || cfg->method == rtBLBadouel || cfg->method == rtBLBadouelGrid) { //Plucker ray-tracing
        pnorm.x = (&(tracer->n[offs].x))[faceid];
        pnorm.y = (&(tracer->n[offs].x))[faceid + 4];
        pnorm.z = (&(tracer->n[offs].x))[faceid + 8];
    } else if (cfg->method == rtHavel || cfg->method == rtBadouel) {
        pn = tracer->m + (offs + faceid) * 3;
    }

    /*pn pointing outward*/

    /*compute the cos of the incidence angle*/
    Icos = fabs(vec_dot(c0, pn));

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
            *eid = *oldeid; /*stay with the current element*/

            if (cfg->debuglevel & dlReflect) {
                MMC_FPRINTF(cfg->flog, "{\"implicit-reflect\": {\"norm\":[%f, %f, %f], \"v\":[%f, %f, %f], \"eid\":%d, \"inroi\":%d, \"faceid\":%d}}\n", pn->x, pn->y, pn->z, c0->x, c0->y, c0->z, *eid, inroi, faceid);
            }
        } else if (cfg->isspecular == 2 && *eid == 0) {
            // if do transmission, but next neighbor is 0, terminate
        } else {                             /*do transmission*/
            vec_mult_add(pn, c0, -Icos, 1.f, c0);
            vec_mult_add(pn, c0, tmp2, n1 / n2, c0);

            if (cfg->debuglevel & dlReflect) {
                MMC_FPRINTF(cfg->flog, "{\"implicit-transmit\": {\"norm\":[%f, %f, %f], \"v\":[%f, %f, %f], \"eid\":%d, \"inroi\":%d, \"faceid\":%d}}\n", pn->x, pn->y, pn->z, c0->x, c0->y, c0->z, *eid, inroi, faceid);
            }
        }
    } else { /*total internal reflection*/
        vec_mult_add(pn, c0, -2.f * Icos, 1.f, c0);
        *eid = *oldeid;

        if (cfg->debuglevel & dlReflect) {
            MMC_FPRINTF(cfg->flog, "{\"implicit-totalreflect\": {\"norm\":[%f, %f, %f], \"v\":[%f, %f, %f], \"eid\":%d, \"inroi\":%d, \"faceid\":%d}}\n", pn->x, pn->y, pn->z, c0->x, c0->y, c0->z, *eid, inroi, faceid);
        }
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
    float3 origin = {r->p0.x, r->p0.y, r->p0.z};

    r->slen = rand_next_scatlen(ran);
    r->inroi = 0;

    if (cfg->srctype == stPencil) { // pencil beam, use the old workflow, except when eid is not given
        if (r->eid > 0) {
            return;
        }
    } else if (cfg->srctype == stPlanar || cfg->srctype == stPattern || cfg->srctype == stPattern3D || cfg->srctype == stFourier || cfg->srctype == stPencilArray) {
        float rx = rand_uniform01(ran);
        float ry = rand_uniform01(ran);

        if (cfg->srctype == stPattern3D) {
            float rz = rand_uniform01(ran);
            r->p0.x = cfg->srcpos.x + rx * cfg->srcparam1.x;
            r->p0.y = cfg->srcpos.y + ry * cfg->srcparam1.y;
            r->p0.z = cfg->srcpos.z + rz * cfg->srcparam1.z;
            r->weight = (cfg->srcnum > 1) ? 1.f :
                        cfg->srcpattern[MIN((int)(rz * JUST_BELOW_ONE * (int)cfg->srcparam1.z), (int)cfg->srcparam1.z - 1) * (int)cfg->srcparam1.y * (int)cfg->srcparam1.x
                                                     + MIN((int)(ry * JUST_BELOW_ONE * (int)cfg->srcparam1.y), (int)cfg->srcparam1.y - 1) * (int)cfg->srcparam1.x
                                                     + MIN((int)(rx * JUST_BELOW_ONE * (int)cfg->srcparam1.x), (int)cfg->srcparam1.x - 1)];
            origin.x += cfg->srcparam1.x * 0.5f;
            origin.y += cfg->srcparam1.y * 0.5f;
            origin.z += cfg->srcparam1.z * 0.5f;
        } else if (cfg->srctype == stPencilArray) {
            r->p0.x = cfg->srcpos.x + floorf(rx * cfg->srcparam1.w) * cfg->srcparam1.x / (cfg->srcparam1.w - 1.f)
                      + floorf(ry * cfg->srcparam2.w) * cfg->srcparam2.x / (cfg->srcparam2.w - 1.f);
            r->p0.y = cfg->srcpos.y + floorf(rx * cfg->srcparam1.w) * cfg->srcparam1.y / (cfg->srcparam1.w - 1.f)
                      + floorf(ry * cfg->srcparam2.w) * cfg->srcparam2.y / (cfg->srcparam2.w - 1.f);
            r->p0.z = cfg->srcpos.z + floorf(rx * cfg->srcparam1.w) * cfg->srcparam1.z / (cfg->srcparam1.w - 1.f)
                      + floorf(ry * cfg->srcparam2.w) * cfg->srcparam2.z / (cfg->srcparam2.w - 1.f);
            r->weight = 1.f;
            origin.x += (cfg->srcparam1.x + cfg->srcparam2.x) * 0.5f;
            origin.y += (cfg->srcparam1.y + cfg->srcparam2.y) * 0.5f;
            origin.z += (cfg->srcparam1.z + cfg->srcparam2.z) * 0.5f;
        } else {
            r->p0.x = cfg->srcpos.x + rx * cfg->srcparam1.x + ry * cfg->srcparam2.x;
            r->p0.y = cfg->srcpos.y + rx * cfg->srcparam1.y + ry * cfg->srcparam2.y;
            r->p0.z = cfg->srcpos.z + rx * cfg->srcparam1.z + ry * cfg->srcparam2.z;
            r->weight = 1.f;

            if (cfg->srctype == stPattern) {
                int xsize = (int)cfg->srcparam1.w;
                int ysize = (int)cfg->srcparam2.w;
                r->posidx = MIN((int)(ry * ysize), ysize - 1) * xsize + MIN((int)(rx * xsize), xsize - 1);
                r->weight = (cfg->srcnum == 1) ? cfg->srcpattern[r->posidx] : 1.f;

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
        }
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
    } else if (cfg->srctype == stDisk || cfg->srctype == stRing || cfg->srctype == stGaussian) { // uniform disk, ring, and Gaussian beam
        float sphi, cphi, phi;

        if (cfg->srctype == stRing && (cfg->srcparam1.z > 0.f || cfg->srcparam1.w > 0.f)) {
            phi = fabsf(cfg->srcparam1.z - cfg->srcparam1.w) * rand_uniform01(ran) + fminf(cfg->srcparam1.z, cfg->srcparam1.w);
        } else {
            phi = TWO_PI * rand_uniform01(ran);
        }

        sphi = sinf(phi);
        cphi = cosf(phi);
        float r0;

        if (cfg->srctype == stDisk || cfg->srctype == stRing) {
            r0 = sqrtf(rand_uniform01(ran) * fabsf(cfg->srcparam1.x * cfg->srcparam1.x - cfg->srcparam1.y * cfg->srcparam1.y) + cfg->srcparam1.y * cfg->srcparam1.y);
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

        if (cfg->srctype == stCone) { // a solid-angle section of a uniform sphere, sampled directly (a rejection loop never ended for a zero angle)
            ang = (cfg->srcparam1.y > 0) ? fminf(fmaxf(cfg->srcparam1.x, 0.f), TWO_PI) * rand_uniform01(ran) :   // uniform in angle
                  acosf(1.f - rand_uniform01(ran) * (1.f - cosf(fminf(fmaxf(cfg->srcparam1.x, 0.f), M_PI))));   // uniform in solid angle
        } else {
            if (cfg->srctype == stIsotropic) { // uniform sphere
                ang = acosf(2.f * rand_uniform01(ran) - 1.f);    //sine distribution
            } else {
                ang = M_PI * rand_uniform01(ran);    //uniform distribution in zenith angle, arcsine
            }
        }

        stheta = sinf(ang);
        ctheta = cosf(ang);
        rotatevector(&(r->vec), stheta, ctheta, sphi, cphi);
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
        rotatevector(&(r->vec), stheta, ctheta, sphi, cphi);
        canfocus = 0;
    } else if (cfg->srctype == stHyperboloid) { // hyperboloid Gaussian beam
        float sphi, cphi;
        float r0 = TWO_PI * rand_uniform01(ran); // phi
        sphi = sinf(r0);
        cphi = cosf(r0);
        r0 = sqrtf(0.5f * rand_next_scatlen(ran)) * cfg->srcparam1.x; // beam waist radius

        /** compute direction in beam-local frame */
        float tilt = (cfg->srcparam1.z > 0.f) ? (-cfg->srcparam1.y / cfg->srcparam1.z) : 0.f;
        float bdir_x = -r0 * sphi;
        float bdir_y =  r0 * cphi;
        float bdir_z =  cfg->srcparam1.z;
        float norm_bdir = 1.f / sqrtf(bdir_x * bdir_x + bdir_y * bdir_y + bdir_z * bdir_z);
        bdir_x *= norm_bdir;
        bdir_y *= norm_bdir;
        bdir_z *= norm_bdir;

        float dp_x = r0 * (cphi - tilt * sphi);
        float dp_y = r0 * (sphi + tilt * cphi);

        /** rotate position offset and direction to srcdir frame */
        if (cfg->srcdir.z > -1.f + EPS && cfg->srcdir.z < 1.f - EPS) {
            float sin_theta = sqrtf(1.f - cfg->srcdir.z * cfg->srcdir.z);
            float inv_sin_theta = 1.f / sin_theta;
            float cphi_beam = cfg->srcdir.x * inv_sin_theta;
            float sphi_beam = cfg->srcdir.y * inv_sin_theta;

            r->p0.x = cfg->srcpos.x + dp_x * cphi_beam * cfg->srcdir.z - dp_y * sphi_beam;
            r->p0.y = cfg->srcpos.y + dp_x * sphi_beam * cfg->srcdir.z + dp_y * cphi_beam;
            r->p0.z = cfg->srcpos.z - dp_x * sin_theta;

            r->vec.x = bdir_x * cphi_beam * cfg->srcdir.z - bdir_y * sphi_beam + bdir_z * cphi_beam * sin_theta;
            r->vec.y = bdir_x * sphi_beam * cfg->srcdir.z + bdir_y * cphi_beam + bdir_z * sphi_beam * sin_theta;
            r->vec.z = -bdir_x * sin_theta + bdir_z * cfg->srcdir.z;
        } else {
            r->p0.x = cfg->srcpos.x + dp_x;
            r->p0.y = cfg->srcpos.y + dp_y;
            r->p0.z = cfg->srcpos.z;

            r->vec.x = bdir_x;
            r->vec.y = bdir_y;
            r->vec.z = (cfg->srcdir.z > 0.f) ? bdir_z : -bdir_z;
        }

        float norm_vec = 1.f / sqrtf(r->vec.x * r->vec.x + r->vec.y * r->vec.y + r->vec.z * r->vec.z);
        r->vec.x *= norm_vec;
        r->vec.y *= norm_vec;
        r->vec.z *= norm_vec;
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
            float4 vv;
            vv.x = r->vec.y * p - r->vec.z * s;
            vv.y = r->vec.z * t - r->vec.x * p;
            vv.z = r->vec.x * s - r->vec.y * t;
            r->vec = vv;
            //*((float3*)&(r->vec))=(float3)(r->vec.y*p-r->vec.z*s,r->vec.z*t-r->vec.x*p,r->vec.x*s-r->vec.y*t);
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
    FLOAT3 vecS = {0.f}, vecAB, vecAC, vecN;
    FLOAT3* nodes = mesh->node;
    int is, i, ea, eb, ec;
    float bary[4] = {0.f};

    for (is = -1; is < mesh->srcelemlen; is++) {
        int include = 1;
        int* elems = NULL;

        if (is < 0) {
            if (r->eid >= 0) {
                elems = (int*)(mesh->elem + (r->eid - 1) * mesh->elemlen);
            } else {
                continue;
            }
        } else {
            elems = (int*)(mesh->elem + (mesh->srcelem[is] - 1) * mesh->elemlen);
        }

        for (i = 0; i < 4; i++) {
            ea = elems[out[i][0]] - 1;
            eb = elems[out[i][1]] - 1;
            ec = elems[out[i][2]] - 1;
            vec_diff3(&nodes[ea], &nodes[eb], &vecAB); //repeated for all photons
            vec_diff3(&nodes[ea], &nodes[ec], &vecAC); //repeated for all photons
            vec_diff3(&nodes[ea], (FLOAT3*) & (r->p0), &vecS);
            vec_cross3(&vecAB, &vecAC, &vecN); //repeated for all photons, vecN can be precomputed
            bary[facemap[i]] = -vec_dot3(&vecS, &vecN);
        }

        for (i = 0; i < 4; i++) {
            if (bary[i] < -1e-4f) {
                include = 0;
            }
        }

        if (include) {
            r->eid = (is >= 0 ? mesh->srcelem[is] : r->eid);
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

    medium* prop = mesh->med + ((cfg->implicit && r->inroi) ? mesh->prop : mesh->type[eid]);
    float ww = r->weight;

    r->weight *= prop->mus / (prop->mua + prop->mus);

    ww -= r->weight;
    r->Eabsorb += ww;

    if (mesh->type[eid] <= 0) { /* no deposit in non-tissue (label<=0) elements */
        return;
    }

    tshift = MIN(MAX((int)((r->photontimer - cfg->tstart) * visit->rtstep), 0), cfg->maxgate - 1) * datalen;

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

/**
 * @brief Retrieve the roisize (radius of node or edge or thickness of face) if present in the current element
 *
 * This function returns ROI size for iMMC simulations
 *
 * \param[in] immctype: cfg->implicit = 1: node or edge iMMC, 2: face iMMC
 * \param[in,out] r: the current ray, r->roisize is the output
 * \param[in] mesh: the mesh data structure
 */

void updateroi(int immctype, ray* r, tetmesh* mesh) {
    if (immctype == 1) {    /* every primitive overlapping the element is in its list: an empty list is outside */
        r->inroi &= (mesh->roioff && mesh->roioff[r->eid] > mesh->roioff[r->eid - 1]);
    } else if (mesh->faceroi) {
        r->roisize = (float*)(mesh->faceroi + (r->eid - 1) * 4);
        // set inroi to 0 if the photon leaves a implicitly marked tetrahedron
        // otherwise do not update
        r->inroi &= (r->roisize[0] != 0.0f);
    }
}
