"""Thin direct subcommands over the pdgstaging operations."""

import argparse
import json

from .operations import stage_source


def _json(text):
    """Parse a CLI JSON argument into Python data.

    Parameters
    ----------
    text : str
        JSON-encoded argument value.

    Returns
    -------
    object
        The decoded Python value.
    """
    return json.loads(text)


def build_parser():
    """Build the CLI parser with one subcommand per operation.

    Returns
    -------
    argparse.ArgumentParser
        Parser with the ``stage-source`` subcommand.
    """
    parser = argparse.ArgumentParser(
        prog="python -m pdgstaging", description="Run a pdgstaging operation"
    )
    subs = parser.add_subparsers(dest="operation", required=True)

    stage = subs.add_parser("stage-source", help="Stage one source into tile shards")
    stage.add_argument("input_path")
    stage.add_argument("shard_root")
    stage.add_argument("source_key")
    stage.add_argument("--z", required=True, type=int, help="Maximum zoom level")
    stage.add_argument("--tms-id", default="WGS1984Quad")
    stage.add_argument(
        "--path-structure",
        default=["style", "tms", "z", "x", "y"],
        type=_json,
        help="JSON list of tile path components",
    )
    stage.add_argument(
        "--properties",
        required=True,
        type=_json,
        help='JSON object of property names, e.g. {"tile": "tile", ...}',
    )
    stage.add_argument("--input-crs", default=None)
    stage.add_argument("--tolerance", default=None, type=float)
    stage.add_argument("--overwrite", action="store_true")

    return parser


def main(argv=None):
    """Parse args, run one operation, and print its result as JSON.

    Parameters
    ----------
    argv : list of str, optional
        Argument list (defaults to ``sys.argv``).

    Returns
    -------
    int
        0 on success; exits with status 2 via ``parser.error`` on
        OSError or ValueError.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = stage_source(
            args.input_path,
            args.shard_root,
            args.source_key,
            tms_id=args.tms_id,
            z=args.z,
            path_structure=tuple(args.path_structure),
            properties=args.properties,
            input_crs=args.input_crs,
            tolerance=args.tolerance,
            overwrite=args.overwrite,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
