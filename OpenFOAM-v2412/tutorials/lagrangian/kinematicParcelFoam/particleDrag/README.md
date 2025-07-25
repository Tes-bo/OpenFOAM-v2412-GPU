Test case for Tomiyama particle drag model

Based on the `parcelInBox` tutorial

- Makes use of the `FreezeParticles` cloud function object to keep the particle
  location fixed in space
- Iterates across a range of particle velocities for different contamination
  levels to reproduce the validation figure from the reference publication

Reference

    Tomiyama, A., Kataoka, I., Zun, I., & Sakaguchi, T. (1998).
    Drag coefficients of single bubbles under normal and micro gravity conditions.
    JSME International Journal Series B Fluids and Thermal Engineering, 41(2), 472-479.
