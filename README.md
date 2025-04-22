# Modelling-the-differences-in-neuronal-activity-caused-by-thermal-vs-millimeter-wave-mediated-heating
We present simulation code modeling neuronal responses to MMW radiation vs. thermal heating. Using a leech interneuron model, we show how active and passive membrane mechanisms influence outcomes. The code enables reproducible studies of non-thermal MMW effects on neurons.


The simulation code was developed to investigate how millimeter wave (MMW) radiation may influence neuronal activity beyond thermal heating effects. Despite ongoing debate in the literature regarding MMW interactions with excitable cells, our goal here is not to resolve the controversy, but to provide a flexible simulation framework that models and compares both thermal and proposed MMW-specific mechanisms affecting neurons.

Using a detailed and physiologically accurate model of a leech ganglion interneuron, we implemented simulations in which the neuron is exposed to a transient linear temperature increase. Traditional thermal effects were governed by the Q10 temperature coefficient, while MMW-specific effects were modeled through dynamic modulation of membrane mechanisms. These included active components such as calcium-dependent potassium channels (IK(Ca)), Na⁺/K⁺-ATPase, and plasma membrane Ca²⁺-ATPase, as well as passive properties like membrane leak conductance, capacitance, and the rate of temperature change (dT/dt).

Our code allows users to vary these parameters systematically across both thermal and MMW scenarios. Simulation results were evaluated based on electrophysiological outputs such as spiking rate and membrane potential dynamics. While thermal heating consistently increased excitability, MMW-modulated simulations often diverged: roughly one-third produced a decrease or complete suppression of spiking activity, reflecting some experimental findings from the literature.

Notably, changes in passive membrane properties had the greatest influence on replicating MMW-specific effects. These results suggest that the impact of MMW radiation may be highly dependent on specific biophysical parameters, reinforcing the importance of experimental model selection.

The codebase we provide supports reproducible, customizable simulations, offering a tool for further investigation into non-thermal biological effects of MMW exposure.

