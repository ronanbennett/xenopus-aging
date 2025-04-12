Example code for loading the clock and predicting on new methylation data:

```python
# Let df_meth_matrix be a pandas DataFrame representing a methylation matrix:
# Rows: CpG sites in the format  Chr:X
# Columns: Samples
# Entries: Methylation value (range from 0 to 1)
from xtropicalis_clock import XTropicalisClock
path_to_clock_weights = "SupplementaryData4.tsv"
clock = XTropicalisClock(path_to_clock_weights)
predicted_ages = clock.predict(df_meth_matrix)
```