---
title: 'MuRAT4: A MATLAB Framework for Multi-Resolution Seismic Attenuation Tomography'
tags:
  - MATLAB
  - seismology
  - attenuation tomography
  - scattering tomography
  - absorption tomography
authors:
  - name: Luca De Siena
    orcid: 0000-0002-3615-5923
    affiliation: "1, 2, 3"
  - name: Chiara Nardoni
    orcid: 0000-0003-1195-3547
    affiliation: "1, 2"
  - name: Martina Avella
    orcid: 0009-0004-2841-5846
    affiliation: "1, 4"
  - name: Gaia Caporale
    orcid: 0009-0006-0261-6839
    affiliation: "4, 5"
  - name: Simona Gabrielli
    orcid: 0000-0002-5829-5900
    affiliation: 6
  - name: Ferdinando Napolitano
    orcid: 0000-0002-2230-078X
    affiliation: 7
  - name: Donato Talone
    orcid: 0009-0004-3472-1022
    affiliation: 5
  - name: Panayiota Sketsiou
    orcid: 0000-0002-6575-4538
    affiliation: 8
  - name: Maria Del Pilar Di Martino
    orcid: 0000-0002-3600-077X
    affiliation: 9
  - name: Thomas King
    orcid: 0000-0002-4536-4073
    affiliation: 10
  - name: Roberto Guardo
    orcid: 0000-0002-2808-9583
    affiliation: "11, 12"
  - name: Edoardo Del Pezzo
    orcid: 0000-0002-6981-5967
    affiliation: "2, 13"
  - name: Francesca Bianco
    orcid: 0000-0001-5400-7724
    affiliation: 2
  - name: Janire Prudencio
    orcid: 0000-0003-1755-3771
    affiliation: 13
  - name: Jesús M. Ibáñez
    orcid: 0000-0002-9846-8781
    affiliation: 13
  - name: Yi Zhang
    orcid: 0000-0002-5421-7737
    affiliation: 14
  - name: Wei-Mou Zhu
    orcid: 0000-0003-4955-543X
    affiliation: 15
  - name: Aqeel Abbas
    orcid: 0000-0002-6133-6857
    affiliation: 16
  - name: Yansong Hu
    orcid: 0000-0002-7512-0678
    affiliation: "1, 17"
  - name: Chiara Saturnino
    orcid: 0009-0009-1124-0556
    affiliation: "1, 4"
  - name: Felix Borleanu
    orcid: 0000-0001-8502-1422
    affiliation: 18
  - name: Liam Bramwell
    orcid: 0009-0005-5475-8882
    affiliation: 3
  - name: Rachit Gautam
    orcid: 0009-0002-9160-6624
    affiliation: "3, 19, 20"
  - name: Miriam C. Reiss
    orcid: 0000-0002-0250-4988
    affiliation: 3
affiliations:
  - name: Dipartimento di Fisica e Astronomia "Augusto Righi", Alma Mater Studiorum, University of Bologna, Bologna, Italy
    index: 1
  - name: Istituto Nazionale di Geofisica e Vulcanologia, Sez. di Napoli - Osservatorio Vesuviano, Napoli, Italy
    index: 2
  - name: Institute of Geosciences, Johannes Gutenberg University Mainz, Mainz, Germany
    index: 3
  - name: Istituto Nazionale di Geofisica e Vulcanologia, Sez. di Bologna, Bologna, Italy
    index: 4
  - name: Dipartimento di Scienze, Università degli studi "G. d'Annunzio", Chieti, Italy
    index: 5
  - name: Istituto Nazionale di Geofisica e Vulcanologia, Sez. di Roma, Roma, Italy
    index: 6
  - name: Università degli Studi di Salerno, Salerno, Italy
    index: 7
  - name: École des Ponts Business School (EPBS), Cyprus
    index: 8
  - name: SLB Milan, Milan, Italy
    index: 9
  - name: Birmingham Institute for Sustainability and Climate Action, University of Birmingham, Birmingham, United Kingdom
    index: 10
  - name: Istituto Nazionale di Geofisica e Vulcanologia, Sez. di Catania - Osservatorio Etneo, Catania, Italy
    index: 11
  - name: Geosciences Barcelona (GEO3BCN), Spanish National Research Council (CSIC), Barcelona, Spain
    index: 12
  - name: Andalusian Institute of Geophysics, University of Granada, Granada, Spain
    index: 13
  - name: Viridien Singapore, Singapore
    index: 14
  - name: Institute of Geology and Geophysics, Chinese Academy of Sciences, Beijing, China
    index: 15
  - name: Department of Earth and Space Sciences, Southern University of Science and Technology, Shenzhen, China
    index: 16
  - name: Institute of Geophysics, China Earthquake Administration, Beijing, China
    index: 17
  - name: National Institute for Earth Physics, Ilfov, Romania
    index: 18
  - name: Ineris ‐ L’Institut national de l’environnement industriel et des risques Campus ARTEM, Nancy Cedex, France
    index: 19
  - name: VITO – Vlaamse Instelling voor Technologisch Onderzoek, Mol, Belgium
    index: 20
date: 4 October 2026
bibliography: paper.bib
---

# Summary

`MuRAT4` (Multi-Resolution Seismic Attenuation Tomography, version 4.0) is an open-access, MATLAB-based framework for measuring and imaging three-dimensional spatial variations of seismic attenuation. It provides tools to separate and image three attenuation mechanisms from proxies measured on seismograms: direct-wave attenuation via the coda-normalization method, scattering attenuation via the peak-delay method, and absorption via late-time coda attenuation. Version 4.0 systematizes these capabilities into a standardized and tested workflow that measures and inverts the three attributes within tomographic frameworks. It is applicable to geological environments ranging from rock volumes to the lithosphere, and has demonstrated applications to seismic-hazard assessment in volcanic, tectonic, geothermal, and subduction settings.

# Statement of Need

Seismic attenuation is a fundamental parameter for understanding subsurface structures and physical processes, constraining heterogeneity, fluid and melt
content, temperature, and seismic energy loss in the Earth's interior, and contributing to seismic-hazard assessment [@mitchell95; @pozgay09; @rawlinson10; @desiena14b]. Attenuation tomography has become an essential complement to velocity imaging for characterizing complex geological environments where scattering and intrinsic absorption strongly modulate seismic amplitudes and waveforms [@desiena09a; @talone23; @napolitano25]. Despite its importance, a critical gap exists in computational seismology for accessible, standardized tools capable of both separating scattering and intrinsic absorption from single-station seismogram amplitudes and imaging attenuation and its mechanisms in a reliable tomographic framework. Existing tools such as Qopen and QEST estimate scattering and intrinsic attenuation at the network scale but do not include direct-wave attenuation as a third observable, and none provides a unified 3D inversion framework for all three parameters. Without such a unified framework, researchers must stitch together incompatible codes to separate the three parameters, then feeding them to alternative inversion frameworks, generally preventing direct comparison of the three parameters. `MuRAT4`, distributed as the [LucaDeSiena/MuRAT](https://github.com/LucaDeSiena/MuRAT) GitHub repository, measures heterogeneous loss from single-station seismograms, separates its mechanisms, and images their spatial variation in 3D using forward tomographic modelling.

# State of the Field

## Open-access codes for measuring and modelling attenuation

Attenuation tomography is established for imaging anelastic structure at global, regional, and local scales [@holiu88; @lees94; @romanowicz1998; @rawlinson10]. However, different mechanisms attenuate seismic waves while they propagate. Several open-source tools are available for estimating or imaging seismic attenuation, although they target different observables and spatial scales. *Qopen* is a Python package that separates intrinsic and scattering attenuation by inversion of seismic-wave envelopes [@eulenfeld16]. Following this research line, the *QEST* code has been developed explicitly for inverting 1-D frequency-dependent scattering and intrinsic attenuation [@vanlaaten24]. *AttenuationVoronoi* implements trans-dimensional Bayesian attenuation tomography using spherical Voronoi cells and has been applied to image attenuation in the Earth's inner core, providing ensembles of models and associated measures of spatial uncertainty [@pejic19].

These packages demonstrate the growing availability of open computational tools for seismic attenuation analysis, while addressing distinct combinations of wave type, attenuation observable, parameterization, and inversion strategy. While existing codes address individual attenuation observables or total attenuation, `MuRAT4` combines multiple seismic attributes to measure attenuation, separate its mechanisms, and model their spatial variations using seismic inversion from earthquake seismograms. MuRAT1.0 [@desiena14a] converted Bash, Fortran, and C++ tools to MATLAB for coda-normalization tomography; MuRAT2.0 added 2D scattering and absorption imaging using peak-delay and Qc proxies [@takahashi07; @calvet13b; @desiena16; @borleanu17; @napolitano20; @guardo22a]; MuRAT3.0 extended the framework to full 3D [@desiena17a; @desiena17b; @reiss22; @napolitano25; @yemane25; @napolitano23; @gabrielli23; @talone23; @talone25; @feng26]. The models shown in Fig. 1 represent results of the application of MuRAT3.0 at Campi Flegrei caldera (Italy) [@desiena17a; @desiena17b]. This and the following 3D figures are obtained by importing the output *.txt* files representing the models into [Paraview](https://www.paraview.org/) using the [GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) [@kaus24].

![High-attenuation (orange) and high-absorption (yellow) isosurfaces at 6 Hz (depth range 0–4 km) from models obtained at Campi Flegrei caldera, Southern Italy, using MuRAT3.0. The high-attenuation isosurface is interpreted as fluid and gas pathways connecting deeper magma to fumaroles; the high-absorption isosurface is interpreted as gas reservoirs feeding the fumaroles.](./Fig1.png){width=100%}

## Related software packages

`MuRAT4` is implemented in MATLAB and uses the Curve Fitting, Signal Processing, Optimization, Global Optimization, Image Processing, and Mapping Toolboxes. Parallel Computing Toolbox is optional. The framework also incorporates open-source routines for SAC file handling, regularization, and geometry and image processing. MATLAB Online access (via the badge above) requires an active institution or individual licence including the listed toolboxes; the full conversion to Julia ([MuRAT.jl](https://github.com/CNardoni/MuRAT.jl.git)) is in development, with current interoperability already available for the GeophysicalModelGenerator.jl visualization step.

# Software Design

`MuRAT4` implements three complementary measurements and several inversion schemes using [SAC files](https://ds.iris.edu/files/sac-manual/) that can be downloaded from open seismological servers such as IRIS/FDSN and must have populated headers:

* **Direct-wave attenuation:** `MuRAT4` uses the coda-normalization (CN) method to estimate direct-wave attenuation from frequency-domain direct and coda-wave energies [@aki80a; @yoshimoto93]. Originally developed for volcanic tomography [@delpezzo06; @desiena08; @desiena09a], its main advantage is in the strong reduction of source and site effects.

* **Scattering attenuation:** Peak delays (PD) quantify the time difference between the direct-wave arrival and maximum envelope amplitude and are used as markers of scattering [@saito02; @takahashi07; @takahashi09; @tripathi10; @calvet13a; @calvet13b; @desiena16; @sato16]. `MuRAT4` accounts for both conventional far-field scattering and strong-scattering effects, including trapped waves and resonances [@napolitano20; @dimartino22b; @desiena16; @gabrielli20; @king22; @king23].

* **Absorption:** `MuRAT4` estimates seismic absorption from coda attenuation (Qc) after a preliminary Multiple-Lapse-Time-Window Analysis (MLTWA) to constrain albedo and extinction lengths [@fehler88; @gabrielli20]. Diffusive modelling and differential sensitivity kernels are used for the inversion [@paasschens97; @pacheco05].

* **Inversion framework:** `MuRAT4` provides standard Tikhonov regularization and alternative iterative MATLAB optimization methods. A regularization matrix $L$ and damping parameter $\epsilon$ control the trade-off between data misfit and model smoothness.

## Example Usage
`MuRAT4` comes with built-in tests and three tutorials as commented live text files that explain the inputs necessary to perform tomography with the three methods and the rationale to select them. In the following, we present the Toba (Indonesia) test case. Two additional test cases are provided for Mount St. Helens (US) and Romania. The Toba test case is based on a dataset of 860 waveforms recorded in the Toba caldera region in 2008. The input file contains the SAC file paths. These files include all source, station and timing information in their header. The inputs are separated by measurement types, inversion, and plotting options (Fig. 2). The output comprises three-dimensional models of direct-wave attenuation, scattering attenuation, and absorption, and several test figures used to assess the results (Fig. 3). The Toba tutorial completes in approximately 130 s on an 8-core MacBook Air M2 (8 GB RAM) running all three methods sequentially; enabling the Parallel Computing Toolbox reduces runtime for the inversion steps.

![Example input file for the application of MuRAT4 to Toba.](./Fig2.png){width=100%}

![Result figures from the Toba test case: the panels show a NW-SE vertical slice through the peak-delay (a), coda attenuation (b), and direct-wave attenuation (c) maps at 3 Hz.](./Fig3.png){width=100%}

# Research Impact Statement

The MuRAT software lineage has been applied across a range of geological environments and spatial scales, including volcanic and geothermal systems, tectonic regions, subduction zones, underground gas storage, and laboratory experiments. The most extensive published record is in volcanic systems. Here, high-attenuation, scattering and absorption anomalies correlate spatially with fluid and gas reservoirs or zones of magma accumulation [@prudencio15a; @prudencio15b; @desiena17a; @desiena17b; @reiss22; @guardo22a]. The relation between high-absorption anomalies and high-temperature gradients and productive geothermal wells has recently opened pathways for applications to geothermal systems [@yemane25; @napolitano25]. At crustal, continental and lithospheric scales, scattering models have also highlighted thrust systems and localized strain near fault tips [@gabrielli23; @talone23], while absorption maps image from sedimentary basins to deeper gas and fluid-saturated volumes [@borleanu23; @talone25; @feng26]. Both peak-delay and coda-attenuation methodologies have been benchmarked and validated in laboratory environments using high-frequency (50 kHz to 1 MHz) acoustic emission and active data from triaxial core experiments [@king23; @dimartino22b].

# AI Usage Disclosure

Claude.ai and MATLAB Copilot were used within MATLAB R2025b to assist with the software development of MuRAT4 . No generative AI tools were used to prepare the documentation or write this manuscript. AI-assisted code was reviewed and validated by the lead author.

# Acknowledgements

This work was funded by Deutsche Forschungsgemeinschaft, Project "A combined imaging and modelling approach to understand magmatic systems across the SE Asia-Australia collision zone", grant 449466527. The lead developer acknowledges the contribution of the original MuRAT contributors, Christine Thomas and Richard Aster.

# References
