/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2024
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
\file    mmc_bench.c

@brief   MMC builtin benchmarks
*******************************************************************************/

#include "mmc_bench.h"

#define MSTR(...) #__VA_ARGS__

const char* benchname[MAX_MCX_BENCH] = {"dmmc-cube60", "dmmc-cube60b", "edgeimmc", "nodeimmc", "faceimmc", ""};

const char* benchjson[MAX_MCX_BENCH] = {
MSTR(
{
    "Session": {
        "ID": "dmmc-cube60",
        "Photons": 1e6,
        "RNGSeed": 1648335518,
        "OutputFormat": "bin",
        "DoMismatch": 0
    },
    "Domain" : {
         "Media": [
             {"mua": 0.00, "mus": 0.0, "g": 1.00, "n": 1.0},
             {"mua": 0.005,"mus": 1.0, "g": 0.01, "n": 1.37}
         ]
    },
    "Mesh": {
        "MeshNode": [
            [ 0, 0, 0],
            [60, 0, 0],
            [ 0,60, 0],
            [60,60, 0],
            [ 0, 0,60],
            [60, 0,60],
            [ 0,60,60],
            [60,60,60]
        ],
        "MeshElem": [
            [1, 2, 8, 4, 1],
            [1, 3, 4, 8, 1],
            [1, 2, 6, 8, 1],
            [1, 5, 8, 6, 1],
            [1, 3, 8, 7, 1],
            [1, 5, 7, 8, 1]
        ],
        "InitElem": 2
    },
    "Forward": {
        "T0": 0.0e+00,
        "T1": 5.0e-09,
        "Dt": 5.0e-09
    },
    "Optode": {
        "Source": {
            "Pos": [30.1, 30.2, 0.0],
            "Dir": [0.0, 0.0, 1.0]
        },
        "Detector": [
            {
                "Pos": [30.0, 20.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [30.0, 40.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [20.0, 30.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [40.0, 30.0, 0.0],
                "R": 1.0
            }
        ]
    }
}),

MSTR(
{
    "Session": {
        "ID": "dmmc-cube60b",
        "Photons": 1e6,
        "RNGSeed": 1648335518,
        "OutputFormat": "bin",
        "DoMismatch": 1
    },
    "Domain" : {
         "Media": [
             {"mua": 0.00, "mus": 0.0, "g": 1.00, "n": 1.0},
             {"mua": 0.005,"mus": 1.0, "g": 0.01, "n": 1.37}
         ]
    },
    "Mesh": {
        "MeshNode": [
            [ 0, 0, 0],
            [60, 0, 0],
            [ 0,60, 0],
            [60,60, 0],
            [ 0, 0,60],
            [60, 0,60],
            [ 0,60,60],
            [60,60,60]
        ],
        "MeshElem": [
            [1, 2, 8, 4, 1],
            [1, 3, 4, 8, 1],
            [1, 2, 6, 8, 1],
            [1, 5, 8, 6, 1],
            [1, 3, 8, 7, 1],
            [1, 5, 7, 8, 1]
        ],
        "InitElem": 2
    },
    "Forward": {
        "T0": 0.0e+00,
        "T1": 5.0e-09,
        "Dt": 5.0e-09
    },
    "Optode": {
        "Source": {
            "Pos": [30.1, 30.2, 0.0],
            "Dir": [0.0, 0.0, 1.0]
        },
        "Detector": [
            {
                "Pos": [30.0, 20.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [30.0, 40.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [20.0, 30.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [40.0, 30.0, 0.0],
                "R": 1.0
            }
        ]
    }
}),

MSTR(
{
    "Session":{
        "ID":"edgeimmc",
        "DoMismatch":1,
        "DebugFlag":"TP",
        "RayTracer":"g",
        "OutputFormat": "bin",
        "Photons":1e6
    },
    "Domain":{
        "Step":[0.01,0.01,0.01],
        "Media":[
            {
                "mua":0,
                "mus":0,
                "g":1,
                "n":1
            },
            {
                "mua":0.0458,
                "mus":35.6541,
                "g":0.9,
                "n":1.37
            },
            {
                "mua":23.0543,
                "mus":9.3985,
                "g":0.9,
                "n":1.37
            }
        ]
    },
    "Mesh":{
        "MeshNode":[
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [1,1,0],
            [0,0,1],
            [1,0,1],
            [0,1,1],
            [1,1,1],
            [0.999,0.5,0.5],
            [0.001,0.5,0.5]
        ],
        "MeshElem":[
            [1,3,10,7,1],
            [8,6,2,9,1],
            [8,3,9,4,1],
            [7,1,5,10,1],
            [4,8,2,9,1],
            [6,5,9,8,1],
            [6,1,2,9,1],
            [1,4,2,9,1],
            [3,1,9,4,1],
            [1,3,9,10,1],
            [1,10,9,5,1],
            [1,5,9,6,1],
            [7,9,10,5,1],
            [9,7,8,5,1],
            [7,9,8,3,1],
            [7,9,3,10,1]
        ],
        "MeshROI":[
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0.1],
            [0,0,0,0.1,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0.1,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0.1,0]
        ],
        "InitElem":14
    },
    "Forward":{
        "T0":0,
        "T1":5e-09,
        "Dt":5e-09
    },
    "Optode":{
        "Source":{
            "Pos":[0.5,0.5,1],
            "Dir":[0,0,-1]
        }
    }
}),

MSTR(
{
     "Session":{
         "ID":"nodeimmc",
         "DoMismatch":1,
         "DebugFlag":"TP",
         "RayTracer":"g",
        "OutputFormat": "bin",
         "Photons":1e6
     },
     "Domain":{
        "Step":[0.01,0.01,0.01],
         "Media":[
             {
                 "mua":0,
                 "mus":0,
                 "g":1,
                 "n":1
             },
             {
                 "mua":0.0458,
                 "mus":35.6541,
                 "g":0.9,
                 "n":1.37
             },
             {
                 "mua":23.0543,
                 "mus":9.3985,
                 "g":0.9,
                 "n":1.37
             }
         ]
     },
     "Mesh":{
         "MeshNode":[
             [0,0,0],
             [1,0,0],
             [0,1,0],
             [1,1,0],
             [0,0,1],
             [1,0,1],
             [0,1,1],
             [1,1,1],
             [0.5,0.5,0.5]
         ],
         "MeshElem":[
             [1,6,9,2,1],
             [4,1,9,2,1],
             [8,4,9,2,1],
             [6,8,9,2,1],
             [6,1,9,5,1],
             [6,5,9,8,1],
             [1,7,9,5,1],
             [5,7,9,8,1],
             [1,3,9,7,1],
             [7,3,9,8,1],
             [3,1,9,4,1],
             [8,3,9,4,1]
         ],
         "MeshROI":[
             [0],
             [0],
             [0],
             [0],
             [0],
             [0],
             [0],
             [0],
             [0.1]
         ],
         "InitElem":8
     },
     "Forward":{
         "T0":0,
         "T1":5e-09,
         "Dt":5e-09
     },
     "Optode":{
         "Source":{
             "Pos":[0.5,0.5,1],
             "Dir":[0,0,-1]
         }
     }
}),

MSTR(
{
     "Session":{
         "ID":"faceimmc",
         "DoMismatch":1,
         "DebugFlag":"TP",
         "RayTracer":"g",
        "OutputFormat": "bin",
         "Photons":1e6
     },
     "Domain":{
        "Step":[0.01,0.01,0.01],
         "Media":[
             {
                 "mua":0,
                 "mus":0,
                 "g":1,
                 "n":1
             },
             {
                 "mua":0.0458,
                 "mus":35.6541,
                 "g":0.9,
                 "n":1.37
             },
             {
                 "mua":23.0543,
                 "mus":9.3985,
                 "g":0.9,
                 "n":1.37
             }
         ]
     },
     "Mesh":{
         "MeshNode":[
             [0,0,0],
             [1,0,0],
             [0,1,0],
             [1,1,0],
             [0,0,1],
             [1,0,1],
             [0,1,1],
             [1,1,1]
         ],
         "MeshElem":[
             [1,2,8,4,1],
             [1,3,4,8,1],
             [1,2,6,8,1],
             [1,5,8,6,1],
             [1,3,8,7,1],
             [1,5,7,8,1]
         ],
         "MeshROI":[
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0.1],
             [0,0,0,0],
             [0,0,0,0.1]
         ],
         "InitElem":6
     },
     "Forward":{
         "T0":0,
         "T1":5e-09,
         "Dt":5e-09
     },
     "Optode":{
         "Source":{
             "Pos":[0.5,0.5,1],
             "Dir":[0,0,-1]
         }
     }
})

};