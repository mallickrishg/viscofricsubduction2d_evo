# Subduction zone earthquake cycles/sequences
 
Impose periodic/aperiodic earthquakes within the coseismic domain and see how the system around the coseismic domain (frictional+viscous domain outside) responds to this stress perturbation by accelerated/decelerated creep. The geometry is a two-dimensional subduction zone setup (plane strain), and the deforming areas of the mantle have been collapsed onto a zero-width shear zone (a fault) to ease computations.

The response of a fault to stress change/evolution is described by rate-strengthening friction (shear_stress = normal_stress*{f_0 + (a-b)*log(v/v_0)} ) while the behaviour of the oceanic and arc mantle are described by viscous rheology. The viscous rheology is defined through the following power-law relationship: strain-rate = A*(stress)^n. Non-local stress interactions are facilitated by a kernel representation of linear elasticity i.e., the boundary element/integral method.

