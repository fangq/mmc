# pmmc examples

Python ports of the MATLAB/Octave demo scripts in `mmclab/example`. Each `demo_*.py`
(and `test_albedo.py`) mirrors the `.m` script of the same name: the same mesh, settings,
printed summaries and roughly the same figures.

All simulations go through `pmmc.mmclab(cfg)`, which applies the same pre- and post-processing
as `mmclab.m`, while `pmmc.run(cfg)` only calls the binary module. Before a simulation it:

- reorients the elements, and computes the face neighbors and element volumes if needed
- finds the initial element, raytracing to the mesh surface when the source is outside the mesh
- adds a source domain to the mesh for wide-field sources (`mmcsrcdomain` and `mmcaddsrc`)

After the simulation it splits the detected photons into fields, sorts the trajectories, and
returns RF outputs as complex arrays. `res['cfg']` holds the prepared mesh; nodal outputs
follow it, as the 3rd output of `mmclab.m` does.

`mmclab.m` saves detected photons, seeds and trajectories depending on how many outputs are
requested. The ports set `cfg['issavedet']`, `cfg['issaveseed']` and `cfg['debuglevel']`
(with `'M'`) explicitly to match.

## Running

```
cd pmmc/example
python3 demo_mmclab_basic.py
```

`genT5mesh`/`genT6mesh` from `mmc/matlab` are replaced by `iso2mesh.meshgrid5`/`meshgrid6`,
which give the same nodes and tetrahedra. Some demos load data with relative paths, so run
them from this folder of an mmc source tree:

- `demo_head_atlas`, `demo_immc_vessel`: `head_atlas.mat` and `vessel.mat` in `mmclab/example`
- `demo_mmcl_b3`: `MMC_Collins_Atlas_Mesh_Version_2L.mat` in `mmclab/example`, downloaded
  separately (see `demo_mmcl_b3.m`)
- `demo_example_meshtest`: `examples/meshtest/sphdiffsemiinf.mat`

## Dependencies

- all demos: `pmmc`, `numpy`, `matplotlib`, `iso2mesh` 0.6.7 or newer (pip install iso2mesh);
  older versions build meshes that differ from MATLAB's or fail in tetgen, and lack the
  `dovol` option of `meshasphere`
- `scipy`: `.mat` loading and interpolation
- `pmcx`: `demo_compare_mmc_mcx`, `demo_dmmc_sphshells`, `demo_dualmesh_output`
- `pmcxcl`: `demo_mmcl_b2_b2d`

`meshplot.py` holds the plane-cut plotting helper shared by several demos.

## Notes

- `demo_optix_mmc` sets `cfg['compute']='optix'`. As in mmclab, a build without OptiX runs
  that case with OpenCL.
- `demo_mmcl_b4` does not run in MATLAB or Octave either. In Octave, as in Python, tetgen
  reports an invalid PLC when `mmcaddsrc` adds the disk source domain. In MATLAB, that step
  succeeds, but mmclab then rejects the mesh for degenerate elements.

## Validation

Every `pmmc.mmclab` call was checked against the matching `mmclab` call. The references come
from Octave, except for `demo_immc_vessel` (which uses Java) and `demo_sfdi_2layer`, which come
from MATLAB. Each pair used the same configuration with fewer photons, and the flux and output
sizes, diffuse reflectance, Jacobian sums and detected photon counts were compared within
Monte Carlo noise.
