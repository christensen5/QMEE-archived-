#!/user/bin/bash

# This script runs and profiles both LV1.py and LV2.py to aid in optimisation.

# Profile the functions.
python -m cProfile LV1.py &

python -m cProfile LV2.py 1. 0.1 1.5 0.75 35 &

python -m cProfile LV3.py &

python -m cProfile LV4.py 1. 0.1 1.5 0.8 35 &

python -m cProfile LV5.py 1. 0.1 1.5 0.8 35 &
