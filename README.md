Many years ago (back in the 90s) I started a simulation program to try out several bone remodeling ideas.
The idea was to use density-dependent material properties in cancellous (spongy) bone - the bone near joint surfces, 
and to adjust that density to mimic the process of bone becoming more dense at areas of high stress,and less dense at areas of low stress.

This is the standard "Wolff's law" hypotheses.

The (continuum-level) mechanical factor that stimulates bone density changes at the time was 
assumed to be the strain energy density divided by the volume fraction of calcified material withing the spongy bone tissue.
Physically this corresponds to using the average strain energy density in the solid phase of the bone tissue.

Simulations adjusting the bone density over time were run by several groups, where density is increased or decreased, 
leading to changes in bone elastic modulus, and that changed modulus lead to changes in the local stress state as the simulations progressed.

These adaptive time-stepping simulations used finite elment models, and the simulations lead to density distributions that approached a continuous distribution, 
but that diverged to states where adjacent individual finite elements were either fully dense or fully empty - the density distribution "checkerboarded".

Since the resulting "checkerboarded" simulations depended on the pattern of the finite elements used in the models, this reslult had some conceptual problems.

One hypothetical way to get a simulation that behaved well (got a contiuous distribution) was to assume that the bone did not react to a 
stimulus when the stimlus was within a range of an equilibrium state - there was a "lazy zone" where remodeling did not react.

Another hypothetical way to arrive at a simulation turned out to be the assumption that the response to stress was more complicated, 
yet still descriptive of "Wolff's law" where higher stress states lead to higher bone density.  If the mechanical stimulus was assumed to be 
the continuum-level strain energy density divided by the calcified mass fraction taken to a power (think mass fraction cubed or to some non-integer power)
then the stimulus would be zero at equilibrium, as in the previous hypothesis, but the slope of the stimulus curve at equilibrium 
could be made zero or close to it. This would be similar to a "lazy zone" but described with a continuous curve.

With that continuous power-law curve describing the remodeling stimulus, 
and with another power-law relationship between bone density and elastic modulus (assuming an isotropic material), some analyis could be done.
The stiffness matrix in the finite elemnt model can be expressed as a sum of element-level matrices (for a density of 1) 
with each element matrix multiplied by the element density taken to a power (see the reference in the repository).

Two things came out - 
1) if the exponent in the stimulus curve is larger than the exponent on the modulus curve, 
the simulations give a continuous, stable, and unique density distribution, for a given set of input parameters.
2) The simulation optimizes a weighted sum of the total strain energy and teh bone density taken to a power.

This lets you do some analytical predictions of the simulation behavior that can be used in prosthesis design.  See the attached paper.

The code and a sample input are still in the process of being fully checked out - they have been rehabilitated after being dormant for a long time.

