# Fresnel Diffraction Patterns

> Numerical computation of near- and far-field Fresnel diffraction patterns through 1D, rectangular, and circular apertures, with a Monte Carlo cross-check.

## Overview

Direct numerical evaluation of the Fresnel diffraction integral for several aperture geometries. Provides a deterministic SciPy `dblquad` solution and an independent Monte Carlo estimator for the circular case, which serves as a cross-check that the two methods converge on the same physical answer.

## Key Features

- **Four diffraction scenarios** in a single interactive menu — 1D slit (near-field), 2D rectangular aperture (far-field), 2D circular aperture (far-field), 2D circular aperture via Monte Carlo
- **Two independent integration strategies** — adaptive deterministic quadrature and stochastic Monte Carlo — for cross-validation
- **Variable y-integration limits** as functions of x to describe the circular aperture boundary in `dblquad`
- **Physically scaled intensity** — output is in units of irradiance (W/m²) via the standard ε₀·c factor on the squared field

## Tech Stack

`Python` · `NumPy` · `SciPy` (`integrate.dblquad`) · `Matplotlib`

## Approach

The Fresnel integral is split into its real and imaginary parts, which are integrated separately and then recombined to give the complex field amplitude. Intensity follows as |E|² scaled by ε₀·c.

For the rectangular case the integration limits are constants and SciPy's `dblquad` handles it directly. The circular case is more interesting: the y-integration limits depend on x via the circle equation, so the limits are passed in as callable functions rather than scalars. This is a small but neat feature of `dblquad` that not every implementation gets right.

The Monte Carlo estimator is included as a sanity check on the deterministic result. It samples random points uniformly inside the aperture's bounding box, rejects those outside the circle, evaluates the Fresnel integrand at each, and averages. Multiplying by the aperture area gives the integral estimate. With 2000 samples per screen point the noise is visible but the central diffraction features are clearly reproduced — the two methods agree to within Monte Carlo statistical error, which is the validation we want.

## Results

The four menu options produce:
- 1D intensity vs screen position for a slit aperture
- 2D intensity map for a rectangular aperture (classic sinc² × sinc² pattern)
- 2D intensity map for a circular aperture (Airy-like rings)
- 2D Monte Carlo estimate of the same circular pattern for comparison

[Add a screenshot of the circular diffraction pattern here — the Airy rings are visually striking.]

## How to Run

```bash
git clone https://github.com/<your-username>/fresnel-diffraction-patterns.git
cd fresnel-diffraction-patterns
pip install -r requirements.txt
python fresnel_diffraction_patterns.py
```

You'll be prompted to choose a problem (`1`, `2`, `3`, `4`, or `q` to quit). Wavelength, screen distance, and aperture sizes are defined as module-level variables and can be changed in-place.

No external data required.
