# ![BDisposalLogos](assets/BDisposal_logo_30x30.png) BDisposal.jl Documentation

Welcome to the documentation for the [_Non parametric productivity analysis under the B-disposal assumption_](https://github.com/sylvaticus/BDisposal.jl) package.

The BDisposal package proposes a serie of environmental efficiency and productivity algorithms for non-parametric modelling when we relax the disposability assumption of some of the outputs and/or inputs (e.g. pollution). These efficiency and productivity measures are implemented through convex and non convex [Data Envelopment Analysis (DEA)](https://en.wikipedia.org/wiki/Data_envelopment_analysis) (aka _Frontier Efficiency Analysis_) models.


## Installation

Until the BDisposal package is registered, please use:
* `] add git@github.com/sylvaticus/BDisposal.jl.git`
Once the BDisposal package will be included in the standard Julia register, install it with:
* `] add BDisposal.jl`

## Loading the module(s)

```julia
using BDisposal
```

## Usage

The BDisposal package contains two prominent components.
First, the BDisposal package defines environmental efficiency indicators for a set of Decision Making Units (`efficiencyScores()`), computing the distance to the best environmental production procedures, i.e., to the efficient production frontier.
Second, the BDisposal package displays environmental productivity indices (`prodIndex()`). These productivity measures are implemented for different time periods (eg., years, months etc.) or spatial units (eg., countries, cities etc.), based on the aforementioned environmental efficiency indicators.

Both components can consider constant or variable returns to scale, convex or non-convex production frontiers and additive or multiplicative distance functions, and work with multiple "good" ("desirable") and "bad" ("undesirable") inputs and outputs (with undesirable inputs being optional).

They are detailed in their respective pages:

- [**`efficiencyScores()`**](efficiencyScores.html): Compute efficiency indicators and convexity test results;
- [**`prodIndex()`**](prodIndex.html): Compute productivity indexes;

## Examples

```julia
# Load Modules
using BDisposal
[....TODO...]
```
## Other packages

This is a list of other Julia packages that use classical DEA method, although without relaxing the disposability assumption.

- [DataEnvelopmentAnalysis.jl](https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl)
- [FrontierEfficiencyAnalysis.jl](https://github.com/wen-chih/FrontierEfficiencyAnalysis.jl)
- [SearchRef.jl](https://github.com/wen-chih/SearchRef.jl)
- [JuMP4DEA.jl](https://github.com/henry8527/JuMP4DEA.jl)


## References

- Abad, A. (2020) [Environmental Efficiency and Productivity Analysis](https://hal.inrae.fr/hal-03032038), _HAL_, 03032038
- Abad, A., P. Ravelojaona (2020a) [Pollution-adjusted Productivity Analysis: The Use of Malmquist and Luenberger Productivity Measures](https://doi.org/10.1002/mde.3260), _Managerial and Decision Economics_
- Abad, A., P. Ravelojoana (2020b) [A Generalization of Environmental Productivity Analysis](https://hal.inrae.fr/hal-02964799), _HAL_, 02964799
- Abad, A., W. Briec (2019) [On the Axiomatic of Pollution-generating Technologies: a Non-Parametric Approach](https://doi.org/10.1016/j.ejor.2019.02.027), _European Journal of Operational Research_, **277**(1), 377-390
- Abad, A. (2018) [Les Enseignements de la Micro-économie de la Production face aux Enjeux Environnementaux: Etude des Productions Jointes. Théorie et Applications](https://tel.archives-ouvertes.fr/tel-01963415), _Ph.D dissertation_, University of Perpignan.
- Abad, A., P. Ravelojaona (2017) [Exponential environmental productivity index and indicators](https://doi.org/10.1007/s11123-017-0513-7), _Journal of Productivity Analysis_, **48**(2), 147-166.
- Abad, A. (2015) [An environmental generalised Luenberger-Hicks-Moorsteen productivity indicator and an environmental generalised Hicks-Moorsteen productivity index](https://doi.org/10.1016/j.jenvman.2015.06.055), _Journal of Environmental Management_, **161**, 325-334.






## Acknowledgements

The development of this package at the _Bureau d'Economie Théorique et Appliquée_ (BETA, Nancy) was supported by the French National Research Agency through the [Laboratory of Excellence ARBRE](http://mycor.nancy.inra.fr/ARBRE/), a part of the “Investissements d'Avenir” Program (ANR 11 – LABX-0002-01).

[![BLogos](assets/logos_betaumr.png)](http://www.beta-umr7522.fr/)
