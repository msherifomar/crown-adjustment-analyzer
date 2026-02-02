# Crown Adjustment Analyzer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A Python tool for quantifying chairside adjustments in dental crown prosthetics using 3D mesh analysis.

## Overview

Crown Adjustment Analyzer compares pre-adjustment (baseline) and post-adjustment STL files to quantify material removed during chairside crown adjustments. Outputs include volume calculations, surface deviation heat maps, and statistical metrics (RMS, mean, SD) for research publication.

## Installation

```bash
git clone https://github.com/mohamedsherifo/crown-adjustment-analyzer.git
cd crown-adjustment-analyzer
pip install -r requirements.txt
python crown_analyzer.py
```

## Usage

1. Launch: `python crown_analyzer.py`
2. Load baseline STL (pre-adjustment)
3. Load adjusted STL (post-adjustment)  
4. Click "Run Analysis"
5. Export results as CSV, report, or screenshot

### Batch Processing

For multiple crown pairs:
```bash
python batch_analyze.py --manifest pairs.csv --output results.csv
```

## Output Metrics

| Metric | Unit |
|--------|------|
| Delta Volume | mm³ |
| Percent Volume Removed | % |
| RMS Deviation | µm |
| Mean Deviation | µm |
| Adjustment Area | mm² |

## Validation

Validated against MeshLab and established mesh analysis methods using 28 crown pairs:
- Volume: Pearson r = 1.000 (p < 0.001)
- RMS deviation: Pearson r = 0.997 (p < 0.001)

## Citation

```bibtex
@article{omar2025crown,
  title={Crown Adjustment Analyzer: A Python Tool for Quantifying Chairside Adjustments in Dental Crown Prosthetics},
  author={Omar, Mohamed Sherif},
  journal={Journal of Open Source Software},
  year={2025},
  doi={10.21105/joss.XXXXX}
}
```

## License

MIT License

## Author

Mohamed Sherif Omar, BDS, MSD  
Clinical Assistant Professor of Prosthodontics  
Director, Digital Innovation Laboratory  
Indiana University School of Dentistry
