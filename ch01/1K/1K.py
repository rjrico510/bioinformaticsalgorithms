# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1K
   Compute the freqeucncy array of a string

"""
import argparse


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Given a text string and k, produce a frequency array"
    )
    parser.add_argument("data_file", help="2-line file - 1st is txt, 2nd is k")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line is txt, 2nd is k

    Returns:
        tuple: (txt, k)
    """
    with open(filename) as f:
        lines = f.readlines()
        txt = lines[0].strip()
        k = int(lines[1].strip())
    return (txt, k)


def symbol_to_number(c: str) -> int:
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


def frequency_array(txt: str, k: int) -> list:
    """compute the frequency array of a text string for strings of length k (k-mers)

    Args:
        text (str): text
        k (int): k-mer length

    Returns:
        list: frequency array
    """
    result = [0] * 4**k

    len_txt = len(txt)
    for i in range(0, len_txt - k + 1):
        result[pattern_to_number(txt[i : i + k])] += 1

    return result


def main():
    """main"""
    args = parse_arguments()
    (txt, k) = parse_file(args.data_file)

    result = frequency_array(txt, k)
    print(f"{' '.join(str(i) for i in result)}")


if __name__ == "__main__":
    main()
