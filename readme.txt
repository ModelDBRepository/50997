"A Computational Model of the Ribbon Synapse"
M.A. Sikora, J. Gottesman, R.F. Miller
J Neurosci Methods (2005) 145:47-61
http://dx.doi.org/10.1016/j.jneumeth.2004.11.023

A model of the ribbon synapse was developed to replicate both pre- 
and postsynaptic functions of this glutamatergic juncture. The 
presynaptic portion of the model is rich in anatomical and 
physiological detail and includes multiple release sites for each
ribbon based on anatomical studies of presynaptic terminals, 
presynaptic voltage at the terminal, the activation of voltage-gated
calcium channels and a calcium-dependent release mechanism whose 
rate varies as a function of the calcium concentration that is 
monitored at two different sites which control both  an ultrafast,
docked pool of vesicles and a release ready pool of tethered 
vesicles. The postsynaptic portion of the program models diffusion of
glutamate and the physiological properties of glutamatergic 
neurotransmission in target cells. We demonstrate the behavior of the
model using the retinal bipolar cell to ganglion cell ribbon synapse.
The model was constrained by the anatomy of salamander bipolar
terminals based on the ultrastructure of these synapses and 
presynaptic contacts were placed onto realistic ganglion cell 
morphology activated by a range of ribbon synapses (46-138). 
These inputs could excite the cell in a manner consistent with 
physiological observations. This model is a comprehensive, first 
generation attempt to assemble our present understanding of the 
ribbon synapse into a domain that permits testing our understanding
of this important structure. We believe that with minor 
modifications of this model it can be fine tuned for other ribbons 
synapses.
~

Figures 7-11 of "A Computational Model of the Ribbon Synapse" can be 
replicated with this code. After compiling all MOD files load the 
appropriate HOC file (fig7.hoc, fig8.hoc, or fig9-11.hoc) from the
neuron main menu after starting nrngui or by either double clicking 
in mswin's window explorer or typing nrngui filename.hoc
in unix.

Session files are called by the HOC file, except for fig8, where 
there are 4 different SES files.  The session files contain 
instructions for running the simulation which are displayed in
a message in a graph when run.

Implementation details are best found in ribbon_tiger.mod

20120113 Updated capump.mod casimple.mod from euler to derivimplicit
as per http://www.neuron.yale.edu/phpbb/viewtopic.php?f=28&t=592
