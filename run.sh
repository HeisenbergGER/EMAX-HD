#!/bin/bash
python setup.py build_ext --inplace
python -c "import main"
