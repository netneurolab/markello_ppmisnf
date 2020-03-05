# -*- coding: utf-8 -*-
"""
Contains full filepaths to various directories relevant to the project

Attribute   Description
---------   -----------
project   : project directory (contains all code, scripts, data, etc.)
data      : primary data directory
bids      : raw BIDS-formatted PPMI imaging data
ppmi      : raw PPMI behavioral data
rois      : ROI/parcellation directory
derived   : derived PPMI data directory
ants      : imaging data processed with ANTS
parcels   : parcellated imaging data processed with ANTS
snf       : directory with data processed with (or pre-processed for) SNF
figs      : directory where figures should be saved
"""

import os
import os.path as op

project = op.dirname(op.dirname(op.dirname(op.abspath(__file__))))
data = op.join(project, 'data')
bids = op.join(data, 'raw', 'ppmi', 'bids')
ppmi = op.join(data, 'raw', 'ppmi', 'behavior')
rois = op.join(data, 'raw', 'rois')
derived = op.join(data, 'derivative')
ants = op.join(data, 'derivative', 'antslct')
parcels = op.join(data, 'derivative', 'parcellated')
snf = op.join(data, 'derivative', 'snf')
figs = op.join(project, 'figures')

for fn in [project, data, bids, ppmi, rois, derived, ants, parcels, snf, figs]:
    os.makedirs(fn, exist_ok=True)

del os, op, fn
