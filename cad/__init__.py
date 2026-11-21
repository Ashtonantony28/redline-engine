"""Redline parametric CAD (M1): build123d parts driven by spec/engine.yaml + spec/cad.yaml.

All geometry is built in millimetres (build123d default, so STEP/STL are in mm); every
quantity crossing the module boundary (mass properties) is converted back to SI (ADR 0001).
All numbers are modelled: there is no physical engine.
"""
