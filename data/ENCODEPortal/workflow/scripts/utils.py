#! bin/bash 

import sys, os
import argparse 

def is_unique(entries):
    a = entries.to_numpy()
    return (a[0] == a[1:]).all()

def check_x_and_y(entries):
    columns = ['Biological replicate(s)_x', 'Biological replicate(s)_y']
    return_values = ['File accession_x', 'File accession_y']
    if not is_unique(entries[columns[0]]):
        print(columns)
        return str(return_values[0]), str(return_values[1])
    elif not is_unique(entries[columns[1]]):
        return str(return_values[1]), str(return_values[0])
    return 

