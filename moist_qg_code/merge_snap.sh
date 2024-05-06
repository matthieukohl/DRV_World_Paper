#!/bin/sh
python3 merge.py /snapshots --cleanup
cd snapshots
python3 -m dedalus merge_sets ./snapshots.h5 snapshots_s*.h5 --cleanup
