# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1D
   Start positions of all instanced of a pattern in a string (0-based)

"""
import argparse


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Given a text string and a pattern, find all start positions for the pattern"
    )
    parser.add_argument(
        "data_file", help="2-line file - 1st is pattern, 2nd is the text"
    )
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line is pattern, 2nd is txt

    Returns:
        tuple: (txt, pattern)
    """
    with open(filename) as f:
        lines = f.readlines()
        pattern = lines[0].strip()
        txt = lines[1].strip()
    return (txt, pattern)


def pattern_match(txt: str, pattern: str) -> list:
    """start positions of all instances of pattern in txt

    Args:
        text (str): text to search
        pattern (str): pattern to match

    Returns:
        list: list of starting positions (int)
    """

    starts = []
    pos = 0
    while pos != -1:
        pos = txt.find(pattern, pos)
        if pos != -1:
            starts.append(pos)
            pos += 1

    return starts


def main():
    """main"""
    args = parse_arguments()
    (txt, pattern) = parse_file(args.data_file)
    print(args.data_file, txt, pattern)

    result = pattern_match(txt, pattern)
    print(f"matches: {' '.join(str(e) for e in result)}")


if __name__ == "__main__":
    main()
