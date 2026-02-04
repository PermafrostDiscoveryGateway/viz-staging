import logging
import argparse
from pdgstaging import TileStager
from pdgstaging import H3GridSummaryGenerator

def valid_h3_resolution(value: str) -> int:
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("H3 resolution must be an integer.")
    if ivalue < 1 or ivalue > 15:
        raise argparse.ArgumentTypeError("H3 resolution must be between 1 and 15.")
    return ivalue


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Stages vector files for the PDG viz tiling pipeline"
    )
    p.add_argument(
        "-c",
        "--config",
        help="Path to configuration JSON file",
        default="config.json",
        type=str,
    )

    p.add_argument("--h3-input", help="Input vector file for H3 summary (.shp/.gpkg)")
    p.add_argument("--h3-output", help="Output H3 summary file (.gpkg)")
    p.add_argument("--h3-res", type=valid_h3_resolution, help="H3 resolution (1–15)")
    p.add_argument("--sum-cols", nargs="*", default=[], help="Attribute columns to sum")
    p.add_argument("--mean-cols", nargs="*", default=[], help="Attribute columns to average")
    p.add_argument("--land-polygons", default=None, help="Optional land/coastline polygon dataset")
    p.add_argument("--area-epsg", type=int, default=6933, help="Equal-area EPSG (default: 6933).")

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    stager = TileStager(args.config)
    stager.stage_all()

    if args.h3_output is not None:
        if args.h3_input is None or args.h3_res is None:
            raise SystemExit("If --h3-output is set, you must also set --h3-input and --h3-res.")

        gen = H3GridSummaryGenerator(area_epsg=args.area_epsg, land_polygons_path=args.land_polygons)
        gen.build_h3_summary(
            input_path=args.h3_input,
            output_path=args.h3_output,
            h3_res=args.h3_res,
            attr_to_sum=args.sum_cols,
            attr_to_mean=args.mean_cols,
            land_polygons_path=args.land_polygons,
            area_epsg=args.area_epsg,
        )

    logging.info("Done")
    logging.shutdown()
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
