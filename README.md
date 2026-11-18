# Flying Wing RCS Explorer

A browser tool for seeing where an airframe reflects radar. Sketch a flying-wing planform, or load a STEP model, and it computes monostatic radar cross-section (RCS) against azimuth with physical optics, then shows which edges, surfaces and parts dominate the return.

Everything runs client-side. STEP geometry is never uploaded.

## What it does

- **Parametric planform.** A thin flat plate with adjustable span, sweep, tip chord, tip cut and sawtooth trailing edge. Edges are grouped by direction, and each group's spike azimuths are predicted and marked.
- **STEP model.** The file is tessellated in the browser with OpenCascade (occt-import-js). Every lit facet is integrated exactly and summed coherently across all azimuths, in parallel Web Workers.
- **Hot-spot view.** A 3D render coloured by backscattered power per unit area, averaged over the nose, beam or tail sectors or taken at a single angle, plus a ranking of parts by their own RCS.
- **Comparison.** Freeze any result as a reference and compare sector averages against it.

## Method in one paragraph

The monostatic physical-optics field of a perfectly conducting surface is the integral of `exp(j·2k·k̂·r)` over its lit area, weighted by `n̂·k̂`. For a flat polygon (or a triangular facet) that integral has a closed form as a sum over its edges (Gordon, 1975), so narrow edge spikes are resolved without meshing. RCS is `σ = 4π/λ² · |Σ (n̂·k̂) I|²`. A full write-up of the formulas and their sources is in the accompanying document.

## Validation

Run from `validation/` with Node 18 or later:

| Script | Checks |
| --- | --- |
| `node test_plate.js` | Polygon integral against the analytic rectangle, a brute-force triangle, and the 4πA²/λ² plate at normal incidence |
| `node test_facet.js` | 3D facet integral, phase included, against brute force |
| `node tG.js` | The POFACETS `G(n,w)` helper: as written in Open RCS, with the recursion corrected, and with a stable series |
| `node tcmp.js` | The closed-form edge integral used here against the corrected POFACETS formula over 7 decades of electrical size |

In-browser check: a 10 mm PEC cube gives −38.54 dBsm face-on at 10 GHz against −38.55 dBsm from 4πA²/λ².

## Limits

Physical optics only. There is no edge diffraction, travelling-wave return, multiple bounce or shadowing between parts, so cavities, inlets and dihedral corners read low, and absolute values near grazing are optimistic. All surfaces are treated as perfect conductors. Use it to compare layouts, not to predict a measured number.

## Files

- `index.html`: the whole app (HTML, CSS and JS in one file)
- `three.min.js`, `OrbitControls.js`: three.js r147 (MIT)
- `occt-import-js.js`, `occt-wasm.js`: occt-import-js 0.0.23 with its WebAssembly inlined as base64 (LGPL-2.1, see `third_party_licenses/`)

To run locally, serve the folder over HTTP (for example `python3 -m http.server`) and open `index.html`.

## References

- W. B. Gordon, "Far-field approximations to the Kirchhoff-Helmholtz representations of scattered fields," *IEEE Trans. Antennas Propag.*, vol. 23, no. 4, 1975.
- D. C. Jenn, *Radar and Laser Cross Section Engineering*, 2nd ed., AIAA, 2005; and POFACETS, Naval Postgraduate School.
- F. Chatzigeorgiadis, *Development of Code for a Physical Optics Radar Cross Section Prediction and Analysis Application*, MS thesis, NPS, 2004.
- E. F. Knott, J. F. Shaeffer, M. T. Tuley, *Radar Cross Section*, 2nd ed., SciTech, 2004.
- Open RCS (comp-ime-eb-br/open-rcs), used as a reference implementation of POFACETS.

Built by Shay Colaco.
