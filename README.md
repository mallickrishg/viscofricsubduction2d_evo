# Subduction zone earthquake cycles/sequences
 
Compute the response of a subduction zone to periodic events or sequences imposed at specified times in a kinematically and dynamically consistent manner.

1. The geometry is a two-dimensional subduction zone setup ([plane strain - ESPM](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JB006611)), and the deforming areas of the mantle have been collapsed onto a zero-width shear zone (a fault) to ease computations
2. The stress calculations are done in using the eigen-strain representation of elastic response to equivalent body-forces. These stress interactions are facilitated by a kernel representation $(K)$ of linear elasticity i.e., the boundary element/integral method $\left(\dot{\tau} = K(v - v^{\infty})\right)$
3. Impose periodic/aperiodic earthquakes within the coseismic domain 
4. The domain where stress is computed and evolves is explicitly outside the coseismic domain 
5. The post- and inter-seismic domain is parameterized by rate-dependent friction $\left(\tau = \sigma_n\left(f_0 + (a-b)*\log(\frac{v}{v_0})\right)\right)$ and by power-law viscous flow $\left(v = A\tau^n\right)$ 

To run the code, use the MAIN_simpleSZ_imposedcycles.m file. Change the fault $(a,b,\sigma_n)$ and viscous shear zone $(A^{-1},n)$ properties accordingly - they can be a spatially fixed value or they vary arbitrarily in space. In the latter case, provide them as a vector with the same dimensions as the discretized boudnary elements. Currently the code works for power-laws only for the upper interface in the ESPM, and odd powers for the lower interface. 

For issues, questions and if you want to discuss please write to me rmallick@caltech.edu

More information about this method and possible uses are shown in a recent [AGU 2022 poster](https://www.authorea.com/doi/full/10.1002/essoar.10513000.1) Tobias Köhne, Rishav Mallick, Mark Simons. Investigating the Potential of Multisequence Displacement Timeseries for Fault Rheology Estimation. Authorea. December 05, 2022.
DOI: 10.1002/essoar.10513000.1
