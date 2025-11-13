Usage of the Mosig-Michalski algorithm for computing Sommerfeld integrals.

Given the following integral that considers a current in a circular region of
radius $\tau$, where $z$ and $z'$ are the perpendicular distances to the stratification
for the field and source, respectively, and $\rho$ is the lateral distance between these
points.

$$
Math: \int_0^\infty \tilde{G}(k_\rho; z|z') \, J_\nu(k_\rho \rho) \, J_\sigma(k_\rho \tau) \, k_\rho \, \mathrm{d}k_\rho
$$

These integrals are often oscillatory and divergent, with poles (singularities)
on the real axis at $a_1$, $a_2$, etc. To compute the tail of this integral more
efficiently, partitioning and extrapolation (via the Mosig-Michalski algorithm)
are applied to the infinite series that approximates the integral.

# References
- MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
- MICHALSKI, K. A. On the efficient evaluation of integral arising in the sommerfeld halfspace problem. In: IEE Proceedings H (Microwaves, Antennas and Propagation). IET Digital Library, 1985. p. 312-318.
