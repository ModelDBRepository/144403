A Moth MGC Model-A HH network with quantitative rate reduction
(Buckley. 2011)

We provide the model used in Buckley (2011). It consists of a network
of Hodgkin Huxley neurons coupled by slow GABA_B synapses which is run
alongside a quantitative reduction described in the associated paper.

Buckley CL, Nowotny T (2011) Multiscale model of an inhibitory network
shows optimal properties near bifurcation.
Phys Rev Lett. 106(23):238109

@ARTICLE{buckley11a, author = {Buckley, C. L. and Nowotny, T}, title =
{Multiscale model of an inhibitory network shows optimal properties
near bifurcation}, journal = {Physical Review Letters}, year =
{{2011}}, volume = {106, 238109}, owner = {clb27}, timestamp =
{2011.02.23} }

It was constructed in linux under gcc

I have provided an exclipse project file but is probably best to
cosntruct your own make file to include the libraries. All the
simulation relevant parameters are contained within

MGCNetwork.cc

fullMGC.h 

is the HH implmenetation and sits aside a quantitative rate equivalent in

fullMGCRate.h



Abstract

We present a systematic multi-scale reduction of a biologically
plausible model of the inihibitory neuronal network of the pheromone
system of the moth. Starting from a Hodgkin-Huxley conductance based
model we adiabatically eliminate fast variables and quantitatively
reduce the model to mean field equations. We then prove analytically
that the network’s ability to operate on signal amplitudes across
several orders of magnitude is optimal when a disinhibitory mode is
close to losing stability and the network dynamics are close to
bifurcation. This has the potential to extend the idea that optimal
dynamic range in the brain arises as a critical phenomenon of phase
transitions in excitable media to brain regions that are dominated by
inhibition or have slow dynamics.
