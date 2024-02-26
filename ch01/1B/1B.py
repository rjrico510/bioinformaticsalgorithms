#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1B
   Find the most frequent k-mers in a string

"""
import argparse
import collections

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Given a text string and k, report the most frequent k-mers"
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

def frequency_words(txt: str, k: int) -> tuple:
    """find most frequent instances of strings of length k (k-mers)

    Args:
        text (str): text to search
        k (int): k-mer length

    Returns:
        tuple: most frequent k-mers
    """
    kmers = collections.defaultdict(int)
    len_txt = len(txt)
    for i in range(0,len_txt-k+1):
        kmers[txt[i:i+k]] += 1

    max_count = max(kmers.values())
    result = {key for key, value in kmers.items() if value == max_count}

    return result

def main():
    """main
    """
    args = parse_arguments()
    (txt, k) = parse_file(args.data_file)
    print(args.data_file, txt, k)

    result = frequency_words(txt, k)
    print(f"result: {result}")

if __name__ == "__main__":
    main()