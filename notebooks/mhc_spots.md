---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Parse MHC binding predictions for a set of peptides

## Dependencies

```python
### Tools to be used
import matplotlib
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
```

```python
### Define an antigen sequence

antigen_sequence = "WQERRTASAADFAAAKALAMAMA"
```


```python
### Create 8-mer slices

slices = []

for i in range(8, len(antigen_sequence)):
  slices.append(antigen_sequence[i-8, i])
```