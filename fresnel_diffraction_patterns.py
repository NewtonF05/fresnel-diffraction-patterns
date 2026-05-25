#!/usr/bin/env python
# coding: utf-8

"""
Fresnel Diffraction Patterns from 1D, Rectangular and Circular Apertures
========================================================================

Author: Newton Fernihough

A numerical study of near- and far-field Fresnel diffraction patterns
behind apertures of several shapes, computed by direct numerical
evaluation of the Fresnel diffraction integral.

Four interactive problems are available from a menu:

    1. 1D diffraction through a slit aperture in the near field
    2. 2D diffraction through a rectangular aperture in the far field
    3. 2D diffraction through a circular aperture, with the y-integration
       limits supplied as functions of x to describe the circular boundary
    4. 2D diffraction through a circular aperture computed with a Monte
       Carlo estimator of the aperture integral, included as an
       independent cross-check of the deterministic ``dblquad`` result

The intensity on the screen is obtained from the squared modulus of the
complex Fresnel field, scaled by epsilon_0 * c so that the output has
units of irradiance.

Run the file and pick a problem from the prompt; enter ``q`` to quit.
"""

# Imports
import numpy as np
from scipy.integrate import dblquad  # Used for double integration
import matplotlib.pyplot as plt  # For plotting graphs

"""Functions"""

def fresnel_real(ya, xa, ys, xs, z, k): 
    """Calculates the real part of the Fresnel integral."""
    coeff = (k / (2 * np.pi * z))  # Constant of integration
    phase = (k / (2 * z)) * ((xs - xa)**2 + (ys - ya)**2)  # Cosine calculation
    return coeff * np.cos(phase)  # Real part of the integral

def fresnel_imag(ya, xa, ys, xs, z, k): 
    """Calculates the imaginary part of the Fresnel integral."""
    coeff = (k / (2 * np.pi * z))  #  Constant of integration
    phase = (k / (2 * z)) * ((xs - xa)**2 + (ys - ya)**2)  # Sine calculation
    return coeff * np.sin(phase)  # Imaginary part of the integral

# Defining the y-limits for the circular aperture
def y_lower(x, r): 
    """Lower y-limit of the aperture."""
    return -np.sqrt(r**2 - x**2)

def y_upper(x, r): 
    """Upper y-limit of the aperture."""
    return np.sqrt(r**2 - x**2)

# Monte Carlo simulation for diffraction pattern calculation
def monte_carlo(xs, ys, n, r, k, z): 
    real_sum, imag_sum = 0, 0  # Initialising real and imaginary sums
    for _ in range(n):  # Loop over random samples
        xa = np.random.uniform(-r, r)  # Random x within aperture
        ya = np.random.uniform(-r, r)  # Random y within aperture
        if xa**2 + ya**2 <= r**2:  # Check if (x,y) is inside the aperture
            real_sum += fresnel_real(ya, xa, ys, xs, z, k)
            imag_sum += fresnel_imag(ya, xa, ys, xs, z, k)
    area = np.pi * r**2  # Area of the circular aperture
    return (real_sum / n) * area, (imag_sum / n) * area  # Return average real and imaginary components

"""variables"""

wavelength = 0.5e-6  # Wavelength in m
k = 2 * np.pi / wavelength  # Wave number
c = 3e8  # Speed of light in m/s
z_near = 0.02  # Near field distance
z_far = 0.05  # Far field distance
epsilon_0 = 8.85e-12  # Permittivity of free space

# Part 1 (1D diffraction pattern)
x_vals_1 = np.linspace(-5e-3, 5e-3, 1001)  # x-values for 1D plot
aperture_1 = (3e-5, 5e-5, 3e-5, 5e-5)  # Aperture limits for 1D problem

# Part 2 (Rectangular diffraction pattern)
x_vals_2 = np.linspace(-2e-3, 2e-3, 100)  # x-values for 2D plot
y_vals_2 = np.linspace(-2e-3, 2e-3, 100)  # y-values for 2D plot
aperture_2 = (-2e-5, 2e-5, -2e-5, 2e-5)  # Aperture limits for rectangular aperture
extents_2 = (-1e-2, 1e-2, -1e-2, 1e-2)  # Plot extent for rectangular diffraction

# Part 3 (Circular diffraction pattern)
x_vals_3 = np.linspace(-2e-3, 2e-3, 100)  # x-values for 3D circular diffraction
y_vals_3 = np.linspace(-2e-3, 2e-3, 100)  # y-values for 3D circular diffraction
radius_3 = 2e-5  # Radius for circular aperture
extents_3 = (-2e-3, 2e-3, -2e-3, 2e-3)  # Plot extent for circular diffraction

# Part 4 (Circular diffraction pattern: Monte Carlo)
x_vals_4 = np.linspace(-2e-3, 2e-3, 100)  # x-values for Monte Carlo
y_vals_4 = np.linspace(-2e-3, 2e-3, 100)  # y-values for Monte Carlo
radius_4 = 2e-5  # Radius for circular aperture
n_samples = 2000  # Number of samples for Monte Carlo

"""menu and loops"""

choice = "0"
while choice != "q":  # Main loop to choose problem
    choice = input('Choose a problem: "1", "2", "3", "4" or "q" to quit: ')  # Prompt user for choice
    print("You chose: ", choice)

    if choice == "1":  # Problem 1: 1D diffraction pattern
        print("Problem 1: 1D Diffraction Pattern")
        intensities = []  # List to store intensity values
        for x in x_vals_1:  # Loop over x-values
            real_int, _ = dblquad(fresnel_real, *aperture_1[:2], *aperture_1[2:], args=(0, x, z_near, k))
            imag_int, _ = dblquad(fresnel_imag, *aperture_1[:2], *aperture_1[2:], args=(0, x, z_near, k))
            intensity = epsilon_0 * c * (real_int**2 + imag_int**2)  # Calculate intensity
            intensities.append(intensity)  # Append intensity to the list
        plt.plot(x_vals_1, intensities, label="Result")  # Plot the intensity vs. x
        plt.xlabel("x (m)")
        plt.ylabel("Intensity")
        plt.title("1D Diffraction")
        plt.legend()
        plt.show()

    elif choice == "2":  # Problem 2: Rectangular diffraction pattern
        print("Problem 2: Rectangular Diffraction Pattern")
        intensities = []  # List to store intensity values
        for x in x_vals_2:  # Loop over x-values
            row = []  # Row of intensities for each y
            for y in y_vals_2:  # Loop over y-values
                real_int, _ = dblquad(fresnel_real, *aperture_2[:2], *aperture_2[2:], args=(y, x, z_far, k))
                imag_int, _ = dblquad(fresnel_imag, *aperture_2[:2], *aperture_2[2:], args=(y, x, z_far, k))
                intensity = epsilon_0 * c * (real_int**2 + imag_int**2)  # Calculate intensity
                row.append(intensity)  # Append intensity for this row
            intensities.append(row)  # Add row to the overall intensities list
        plt.imshow(intensities, extent=extents_2, origin="lower", cmap="nipy_spectral_r")  # Display as 2D image
        plt.colorbar(label="Intensity")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title("2D Diffraction")
        plt.show()

    elif choice == "3":  # Problem 3: Circular diffraction pattern
        print("Problem 3: Circular Diffraction Pattern")
        
        # Functions for fixed y-limits based on radius
        def y_lower_fixed(x):
            return y_lower(x, radius_3)

        def y_upper_fixed(x):
            return y_upper(x, radius_3)

        intensities = []  # List to store intensity values
        for x in x_vals_3:  # Loop over x-values
            row = []  # Row of intensities for each y
            for y in y_vals_3:  # Loop over y-values
                real_int, _ = dblquad(fresnel_real, -radius_3, radius_3, y_lower_fixed, y_upper_fixed, args=(y, x, z_far, k))
                imag_int, _ = dblquad(fresnel_imag, -radius_3, radius_3, y_lower_fixed, y_upper_fixed, args=(y, x, z_far, k))
                intensity = epsilon_0 * c * (real_int**2 + imag_int**2)  # Calculate intensity
                row.append(intensity)  # Append intensity for this row
            intensities.append(row)  # Add row to the overall intensities list
        plt.imshow(intensities, extent=extents_3, origin="lower", cmap="nipy_spectral_r")  # Display as 2D image
        plt.colorbar(label="Intensity")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title("2D Diffraction (Circular Aperture)")
        plt.show()

    elif choice == "4":  # Problem 4: Circular Diffraction Pattern: Monte Carlo
        print("Problem 4: Circular Diffraction Pattern: Monte Carlo") 
        intensities = []  # List to store intensity values
        for x in x_vals_4:  # Loop over x-values
            row = []  # Row of intensities for each y
            for y in y_vals_4:  # Loop over y-values
                real_int, imag_int = monte_carlo(x, y, n_samples, radius_4, k, z_far)  # Monte Carlo simulation
                intensity = np.abs(real_int + 1j * imag_int)**2  # Calculate intensity from complex result
                row.append(intensity)  # Append intensity for this row
            intensities.append(row)  # Add row to the overall intensities list
        plt.imshow(intensities, extent=(x_vals_4[0], x_vals_4[-1], y_vals_4[0], y_vals_4[-1]), origin="lower", cmap="nipy_spectral_r")  # Display as 2D image
        plt.colorbar(label="Intensity")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title("Circular Diffraction Pattern (Monte Carlo)")
        plt.show()

    elif choice != "q": 
        print("Invalid choice. Try again.")

print("Goodbye!") 
