---
title: 'Crown Adjustment Analyzer: A Python Tool for Quantifying Chairside Adjustments in Dental Crown Prosthetics'
tags:
  - Python
  - dentistry
  - prosthodontics
  - CAD/CAM
  - 3D mesh analysis
  - dental crowns
authors:
  - name: Mohamed Sherif Omar
    orcid: 0009-0006-9198-3250
    affiliation: 1
affiliations:
  - name: Department of Prosthodontics, Indiana University School of Dentistry, Indianapolis, Indiana, USA
    index: 1
date: 1 February 2025
bibliography: paper.bib
---

# Summary

Chairside adjustment of dental crowns is a critical clinical step that impacts restoration longevity, patient comfort, and procedural efficiency. Despite the emergence of AI-powered crown design platforms promising to reduce clinical adjustment requirements, no standardized open-source tool exists for quantifying the volume and distribution of material removed during these adjustments. Crown Adjustment Analyzer addresses this gap by providing researchers and clinicians with a scientifically validated tool for measuring chairside adjustments using 3D mesh analysis, built upon established open-source libraries and designed for reproducible research workflows.

# Statement of Need

The evaluation of dental crown design platforms—whether human-designed or AI-generated—requires objective metrics for clinical performance. Current literature predominantly focuses on virtual design accuracy (marginal fit, internal adaptation) without addressing the practical clinical outcome of chairside adjustment requirements. This represents a significant gap because: (1) a crown with excellent virtual fit may still require extensive chairside adjustment due to manufacturing tolerances or contact point discrepancies, (2) excessive adjustment can compromise crown integrity and increase fracture risk, and (3) adjustment time directly impacts clinical efficiency and patient experience.

Existing commercial software such as Geomagic Control X provides surface deviation analysis but presents several limitations for research applications: it lacks specific features for dental crown adjustment quantification, requires expensive licenses that limit accessibility, does not provide open and reproducible workflows suitable for research publication, and cannot be extended or validated by the research community.

Crown Adjustment Analyzer addresses these limitations by providing precise volume calculation of material removed during adjustment using the mathematically rigorous signed tetrahedra method, surface deviation heat maps showing adjustment distribution across the crown surface, statistical metrics (RMS, mean, standard deviation) consistent with dental research reporting standards, publication-quality visualizations with exportable data for statistical analysis, and an open-source codebase that enables validation, extension, and reproducibility.

# Design Decisions and Architecture

Crown Adjustment Analyzer was designed with several key architectural decisions to maximize utility for the dental research community.

**Building on established ecosystems.** Rather than implementing mesh processing algorithms from scratch, the software builds upon trimesh [@trimesh], a well-maintained library for 3D mesh operations with over 2,000 GitHub stars and active development. This decision ensures computational reliability, benefits from ongoing community improvements, and allows the software to focus on domain-specific dental analysis features. Visualization leverages PyVista [@sullivan2019pyvista], the standard Python interface to VTK, providing publication-quality 3D rendering.

**Graphical interface for clinical researchers.** While command-line tools offer flexibility for programmers, dental researchers and clinicians often lack programming experience. The PyQt5-based graphical interface provides an intuitive workflow—load files, click analyze, export results—that requires no coding knowledge. This design choice prioritizes accessibility and adoption within the dental research community.

**Standardized output metrics.** The software calculates RMS deviation, which is the standard metric in dental surface comparison studies [@ahrberg2016evaluation; @cho2015evaluating], ensuring that results are directly comparable with existing literature and commercial software benchmarks.

**Extensible architecture.** The core `MeshAnalyzer` class separates analysis logic from the GUI, enabling programmatic use for batch processing, integration into larger research pipelines, or extension with additional analysis methods.

# Methodology

## Volume Calculation

The software employs the signed tetrahedra method for volume calculation, implementing the divergence theorem to compute mesh volume as the sum of signed tetrahedra formed by each triangular face and the coordinate origin [@zhang2001efficient]. For watertight meshes, this provides exact volume measurement; for non-watertight meshes, convex hull approximation is used with appropriate warnings to the user.

## Surface Deviation Analysis

Point-to-surface distances are computed using trimesh's proximity queries, which employ spatial indexing (k-d trees) for efficient nearest-point calculations on meshes with tens of thousands of vertices. The signed distance is determined by:

$$d_{signed} = d_{euclidean} \cdot sign(\vec{v} \cdot \vec{n})$$

Where $\vec{v}$ is the displacement vector from the closest baseline point to the adjusted mesh vertex, and $\vec{n}$ is the interpolated surface normal at the closest point. Negative values indicate material removal (the adjusted surface lies inside the baseline surface), which is the expected direction for chairside adjustments.

## RMS Deviation

Root Mean Square deviation provides a single metric summarizing overall adjustment magnitude:

$$RMS = \sqrt{\frac{1}{n}\sum_{i=1}^{n}d_i^2}$$

This metric is widely used in dental research for surface comparison studies and provides direct comparability with existing literature and commercial software outputs.

# Features

Crown Adjustment Analyzer provides dual-file comparison with automatic mesh validation and repair for non-watertight meshes, volume analysis calculating absolute and percentage volume removed, deviation mapping with color-coded heat maps using multiple colormap options, statistical analysis computing RMS, mean, standard deviation, and percentile statistics, interactive 3D visualization with adjustable lighting, opacity, and view modes, and comprehensive data export including CSV for statistical software, text reports, and publication-ready PNG screenshots.

The graphical interface (\autoref{fig:interface}) provides an intuitive workflow suitable for clinical researchers without programming experience.

![Crown Adjustment Analyzer interface showing heat map visualization of chairside adjustments on a dental crown. The color scale indicates surface deviation in micrometers, with blue regions showing material removal.\label{fig:interface}](figures/interface_screenshot.png)

# Validation

The software was validated using 28 crown pairs across three design groups (Human Expert, 3Shape AI, and Exocad AI). Volume measurements were verified against MeshLab (ISTI-CNR, Italy), with each baseline and adjusted STL file manually opened and volume recorded individually—a total of 56 manual measurements. Crown Adjustment Analyzer volume calculations showed near-perfect agreement with MeshLab (Pearson r = 1.000, p < 0.001, mean absolute difference < 0.01 mm³), confirming the accuracy of the signed tetrahedra implementation. RMS surface deviation measurements were compared against reference values obtained using established mesh analysis methods, showing excellent correlation (Pearson r = 0.997, p < 0.001) with mean absolute difference of 5.3 µm. This validation demonstrates that Crown Adjustment Analyzer produces reliable, reproducible measurements consistent with established open-source tools, making it suitable for quantitative dental research.

# Research Applications

Crown Adjustment Analyzer was developed as part of a comprehensive research program at the Digital Innovation Laboratory at Indiana University School of Dentistry evaluating AI-powered crown design platforms. The software enables standardized comparison of adjustment requirements between human expert and AI-designed crowns, systematic evaluation of different AI platforms (e.g., 3Shape Automate, Dentbird, Medit AI, ExoCAD AI), investigation of correlations between virtual fit metrics and clinical adjustment needs, and assessment of adjustment volume impact on mechanical properties and longevity.

# Related Publications

This software was developed to support an ongoing research study comparing AI-powered dental crown design platforms. A manuscript describing the clinical evaluation—including design time, adjustment time, volume removed, and surface deviation across multiple AI systems—is currently in preparation for submission to the Journal of Prosthodontics. Crown Adjustment Analyzer provides the quantitative analysis methodology for that study, and its development represents a substantial and independently useful contribution to dental research infrastructure.

# Acknowledgements

The author thanks the Digital Innovation Laboratory team at Indiana University School of Dentistry for their support in testing and validating this software.

# References
