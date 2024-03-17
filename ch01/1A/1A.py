# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1A
   Compute the number of times a pattern appears in text

"""
import argparse


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Given a text string and a pattern, count the # instances of the pattern"
    )
    parser.add_argument(
        "data_file", help="2-line file - 1st is txt, 2nd is the pattern"
    )
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line is txt, 2nd is pattern

    Returns:
        tuple: (txt, pattern)
    """
    with open(filename) as f:
        lines = f.readlines()
        txt = lines[0].strip()
        pattern = lines[1].strip()
    return (txt, pattern)


def pattern_count(txt: str, pattern: str) -> int:
    """count instances of pattern in txt

    Args:
        text (str): text to search
        pattern (str): pattern to match

    Returns:
        int: # matches
    """
    count = 0
    len_txt = len(txt)
    len_pattern = len(pattern)
    for i in range(0, len_txt - len_pattern + 1):
        if txt[i : i + len_pattern] == pattern:
            count += 1

    # alt solution using find()
    # pos = 0
    # while pos != -1:
    #     pos = txt.find(pattern, pos)
    #     if pos != -1:
    #         count += 1
    #         pos += 1

    return count


def main():
    """main"""
    args = parse_arguments()
    (txt, pattern) = parse_file(args.data_file)
    print(args.data_file, txt, pattern)

    result = pattern_count(txt, pattern)
    print(f"count: {result}")


if __name__ == "__main__":
    main()
