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

The package provide two functions.

- [**`efficiencyScores()`**](efficiencyScores.html): Compute efficiency scores;
- [**`prodIndex()`**](prodIndex.html): Compute Productivity Indexes;

## Examples

```julia
# Load Modules
using BDisposal
```
## Other packages

This is a list of other Julia packages that use classical DEA method, although without relaxing the disposability assumption.

- [DataEnvelopmentAnalysis.jl](https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl)
- [FrontierEfficiencyAnalysis.jl](https://github.com/wen-chih/FrontierEfficiencyAnalysis.jl)
- [SearchRef.jl](https://github.com/wen-chih/SearchRef.jl)
- [JuMP4DEA.jl](https://github.com/henry8527/JuMP4DEA.jl)

<!-- - [Various code for DEA](https://github.com/wen-chih/code-for-DEA) -->

## References


- Abad, A. (2020) [Environmental Efficiency and Productivity Analysis](https://hal.inrae.fr/hal-03032038), _HAL_, 03032038.
- Abad, A., P. Ravelojaona (2020a) [Pollution-adjusted Productivity Analysis: The Use of Malmquist and Luenberger Productivity Measures](https://doi.org/10.1002/mde.3260), _Managerial and Decision Economics_
- Abad, A., P. Ravelojoana (2020b) A Generalization of Environmental Productivity Analysis, hal-02964799.
- Abad, A., W. Briec (2019) On the Axiomatic of Pollution-generating Technologies: a Non-Parametric Approach, {\it European Journal of Operational Research}, 277(1), 377-390.}\
- Abad, A. (2018) {\it Les Enseignements de la Micro-\'economie de la Production face aux Enjeux Environnementaux: Etude des Productions Jointes. Th\'eorie et Applications}, Ph.D dissertation, University of Perpignan.
- Abad, A., P. Ravelojaona (2017) Exponential environmental productivity index and indicators, {\it Journal of Productivity Analysis}, 48(2), 147-166.\
- Abad, A. (2015) An environmental generalised Luenberger-Hicks-Moorsteen
productivity indicator and an environmental generalised Hicks-Moorsteen productivity index, {\it Journal of Environmental Management}, 161, 325-334.\






## Acknowledgements

The development of this package at the _Bureau d'Economie Théorique et Appliquée_ (BETA, Nancy) was supported by the French National Research Agency through the [Laboratory of Excellence ARBRE](http://mycor.nancy.inra.fr/ARBRE/), a part of the “Investissements d'Avenir” Program (ANR 11 – LABX-0002-01).

[![BLogos](assets/logos_betaumr.png)](http://www.beta-umr7522.fr/)
