"""
MIT License

Console entry-point for TEcounter.
"""

from __future__ import annotations

from .cli import build_parser, dispatch


def main() -> None:
    """Entry-point used by `python -m tecounter` and console script."""
    parser = build_parser()
    args = parser.parse_args()
    dispatch(args)


if __name__ == "__main__":
    main()
