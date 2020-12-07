/***************************************************************************//**
**  \mainpage Monte Carlo eXtreme - GPU accelerated Monte Carlo Photon Migration
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2009-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, 
**          <a href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
**          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated 
**          by Graphics Processing Units,"</a> Optics Express, 17(22) 20178-20190 (2009).
**  \li \c (\b Yu2018) Leiming Yu, Fanny Nina-Paravecino, David Kaeli, and Qianqian Fang,
**          "Scalable and massively parallel Monte Carlo photon transport
**           simulations for heterogeneous computing platforms," J. Biomed. Optics, 23(1), 010504, 2018.
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
