#!/bin/bash
find ./ -name "*.py" -type f -print0 | xargs -0 -I {} static/tidy.py {} {}
