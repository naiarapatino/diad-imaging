#!/usr/bin/env python3
import casatools
print('CASA tools imported successfully!')
tb = casatools.table()
print('Table tool created successfully!')
print(f'CASA version: {casatools.version_string()}')
