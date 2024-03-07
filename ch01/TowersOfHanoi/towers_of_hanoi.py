#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01
   Towers of Hanoi

"""
import argparse
import typing

class peg(typing.NamedTuple):
    start: int
    transit: int
    end: int

PEGS = peg(1,2,3)

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Enumerate moves of Towers of Hanoi for a given # disks (n)"
    )
    parser.add_argument("n", help="# of disks in the tower", type=int)
    args = parser.parse_args()
    return args

def towers_of_hanoi(n: int, start: int, end: int) -> None:
    """Moves for Towers of Hanoi

    Args:
        n (int): # disks
        start (int): starting peg
        end (int): ending peg
    """

    if n == 1:
        print(f"n={n}: move from {start} to {end}")
        return
    
    transit = 6 - start - end
    towers_of_hanoi(n - 1, start, transit)
    print(f"n={n}: move from {start} to {end}")
    towers_of_hanoi(n - 1, transit, end)

def main():
    """main
    """
    args = parse_arguments()

    towers_of_hanoi(args.n, PEGS.start, PEGS.end)

if __name__ == "__main__":
    main()