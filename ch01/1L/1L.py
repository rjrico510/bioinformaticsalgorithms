#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1L
   Implement pattern_to_number

"""
import argparse

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="pattern_to_number"
    )
    parser.add_argument("data_file", help="1-line file of text")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> str:
    """Parse file

    Args:
        filename (str): file - 1-line file of text

    Returns:
        str: txt
    """
    with open(filename) as f:
        lines = f.readlines()
        txt = lines[0].strip()
    return txt

def symbol_to_number(c:str) -> int:
    """convert symbol (A,C,G,T) to int
    
    Args:
        c (str): single character - must be in [ACGT]

    Returns:
        int: index corresponding to pattern
    """
    SYMBOL_MAP = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
    }
    return SYMBOL_MAP[c.upper()]

def pattern_to_number(pattern: str) -> int:
    """Convert pattern to an integer index for the frequency array

    Args:
        pattern (str): pattern

    Returns:
        int: index corresponding to pattern
    """
    if len(pattern) == 0:
        return 0
    
    last_symbol = pattern[-1]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(last_symbol)

def main():
    """main
    """
    args = parse_arguments()
    txt = parse_file(args.data_file)

    result = pattern_to_number(txt)
    print(f"{result}")

if __name__ == "__main__":
    main()