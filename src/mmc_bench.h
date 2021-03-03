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
\file    mcx_bench.h

@brief   MCX builtin benchmarks
*******************************************************************************/
                                                                                         
#ifndef _MCEXTREME_BENCHMARK_H
#define _MCEXTREME_BENCHMARK_H

#define MSTR(...) #__VA_ARGS__

const char *benchname[]={"d-cube60","d-cube60b",""};
const char *benchjson[]={
MSTR(
{
    "Session": {
	"ID":       "d-cube60",
	"Photons":  1e6,
	"RNGSeed":  1648335518,
	"DoMismatch": 0
    },
    "Domain": {
        "Dim":    [60,60,60],
        "OriginType": 1,
        "Media": [
             {"mua": 0, "mus": 0, "g": 1, "n": 1},
             {"mua": 0.005,"mus": 1.0, "g": 0.01, "n": 1.37},
             {"mua": 0.002,"mus": 5, "g": 0.90, "n": 1}
        ]
    },
    "Shapes": [
        {"Name":     "cubic60"},
        {"Origin":   [0,0,0]},
        {"Grid":     {"Tag":1, "Size":[60,60,60]}}
    ],
    "Forward": {
	"T0": 0.0e+00,
	"T1": 5.0e-09,
	"Dt": 5.0e-09
    },
    "Optode": {
	"Source": {
	    "Type":"pencil",
	    "Pos": [29.0, 29.0, 0.0],
	    "Dir": [0.0, 0.0, 1.0]
	},
	"Detector": [
	    {
		"Pos": [29.0,  19.0,  0.0],
		"R": 1.0
	    },
            {
                "Pos": [29.0,  39.0,  0.0],
                "R": 1.0
            },
            {
                "Pos": [19.0,  29.0,  0.0],
                "R": 1.0
            },
            {
                "Pos": [39.0,  29.0,  0.0],
                "R": 1.0
            }
	]
    }
}),


MSTR(
{
    "Session": {
	"ID":       "d-cube60b",
	"Photons":  1e6,
	"RNGSeed":  1648335518,
	"DoMismatch": true
    },
    "Domain": {
        "Dim":    [60,60,60],
        "OriginType": 1,
        "Media": [
             {"mua": 0.00, "mus": 0.0, "g": 1.00, "n": 1.0},
             {"mua": 0.005,"mus": 1.0, "g": 0.01, "n": 1.37},
             {"mua": 0.002,"mus": 5.0, "g": 0.90, "n": 1.0}
        ]
    },
    "Shapes": [
        {"Name":     "cube60b"},
        {"Origin":   [0,0,0]},
        {"Grid":     {"Tag":1, "Size":[60,60,60]}}
    ],
    "Forward": {
	"T0": 0.0e+00,
	"T1": 5.0e-09,
	"Dt": 5.0e-09
    },
    "Optode": {
	"Source": {
	    "Type":"pencil",
	    "Pos": [29.0, 29.0, 0.0],
	    "Dir": [0.0, 0.0, 1.0]
	},
	"Detector": [
	    {
		"Pos": [29.0,  19.0,  0.0],
		"R": 1.0
	    },
            {
                "Pos": [29.0,  39.0,  0.0],
                "R": 1.0
            },
            {
                "Pos": [19.0,  29.0,  0.0],
                "R": 1.0
            },
            {
                "Pos": [39.0,  29.0,  0.0],
                "R": 1.0
            }
	]
    }
})

};

#endif
