#!/usr/bin/env python3
import pandas as pd
import numpy as np
from bootstrap_resample import describe, sample

TEST = dict(
    seed = 42,
    N = int(1e4),
    columns = ['Cats', 'Dogs'],
    metric = np.mean,
    FDR_method = 'Storey'
)

if __name__ == '__main__':
    np.random.seed(seed=TEST['seed'])
    cols = TEST['columns']
    x = pd.DataFrame(np.random.randn(TEST['N'], len(cols)) + np.arange(len(cols)), columns=cols)
    print(describe(x, TEST['metric']))
    print(sample(x, TEST['metric']).pTable(method=TEST['FDR_method']).to_string())
