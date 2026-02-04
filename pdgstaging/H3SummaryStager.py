# H3SummaryStager.py
import logging
import os
from pathlib import Path
from datetime import datetime

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from filelock import FileLock
from typing import Optional

from . import TilePathManager
from .H3GridSummaryGenerator import H3GridSummaryGenerator


class H3SummaryStager:
    """
    Runs H3GridSummaryGenerator over all input vectors and writes one
    H3 summary GeoPackage per input file under a configured base_dir.
    Optionally records a CSV/Parquet run summary.
    """

    def __init__(
        self,
        tiles: Optional[TilePathManager] = None,
        out_base_dir: str = "h3",
        out_ext: str = ".gpkg",
        summary_filename: str = "h3_summary.csv",
        generator: Optional[H3GridSummaryGenerator] = None,
    ):
        self.logger = logging.getLogger(__name__)
        self.tiles = tiles
        self.out_base_dir = out_base_dir
        self.out_ext = out_ext

        if self.tiles is None:
            raise ValueError("tiles (TilePathManager) is required for workflow integration")

        staged_root = self.tiles.base_dirs["staged"]["path"]
        summary_root = os.path.dirname(staged_root)
        self.summary_path = os.path.join(summary_root, summary_filename)

        self.gen = generator or H3GridSummaryGenerator(
            tiles=self.tiles,
            out_base_dir=out_base_dir,
            logger=self.logger,
        )

    def _output_path_for_input(self, input_path: str | Path, h3_res: int) -> Path:
        input_path = Path(input_path)

        input_root = Path(self.tiles.base_dirs["input"]["path"])
        try:
            rel = input_path.relative_to(input_root)
        except Exception:
            rel = Path(input_path.name)

        rel = Path(rel)
        stem = rel.stem
        parent = rel.parent

        # outputs sit next to input/staged directories (same level)
        out_root = Path(self.tiles.base_dirs["input"]["path"]).parent / self.out_base_dir
        out_dir = out_root / parent
        out_dir.mkdir(parents=True, exist_ok=True)

        return out_dir / f"{stem}_h3r{h3_res}{self.out_ext}"

    def stage_all(
        self,
        h3_res: int,
        attr_to_sum=None,
        attr_to_mean=None,
        land_polygons_path: str | Path | None = None,
        area_epsg: int | None = None,
    ):
        overall_start = datetime.now()
        input_paths = self.tiles.get_filenames_from_dir("input")
        n = len(input_paths)

        if n == 0:
            self.logger.error("No vector files found for H3 staging.")
            return

        self.logger.info(f"Begin H3 staging {n} input vector files.")

        rows = []
        for p in input_paths:
            start = datetime.now()
            out = None
            ok = False
            err = None
            try:
                out = self.stage(
                    path=p,
                    h3_res=h3_res,
                    attr_to_sum=attr_to_sum,
                    attr_to_mean=attr_to_mean,
                    land_polygons_path=land_polygons_path,
                    area_epsg=area_epsg,
                )
                ok = True
            except Exception as e:
                err = repr(e)
                self.logger.exception(f"Failed to stage H3 summary for {p}")

            rows.append(
                {
                    "input_path": str(p),
                    "output_path": str(out) if out else None,
                    "h3_res": h3_res,
                    "ok": ok,
                    "seconds": (datetime.now() - start).total_seconds(),
                    "error": err,
                    "datetime": datetime.now().isoformat(),
                }
            )

        df = pd.DataFrame(rows)
        self._append_summary(df)

        total = datetime.now() - overall_start
        self.logger.info(f"H3-staged {n} files in {total} ({total / max(n,1)} per file).")

    def stage(
        self,
        path: str | Path,
        h3_res: int,
        attr_to_sum=None,
        attr_to_mean=None,
        land_polygons_path: str | Path | None = None,
        area_epsg: int | None = None,
        output_path: str | Path | None = None,
    ) -> Path:
        if output_path is None:
            output_path = self._output_path_for_input(path, h3_res)

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        lock = self._lock_file(str(output_path))
        try:
            self.gen.build_h3_summary(
                input_path=path,
                output_path=output_path,
                h3_res=h3_res,
                land_polygons_path=land_polygons_path,
                area_epsg=area_epsg,
                attr_to_sum=attr_to_sum,
                attr_to_mean=attr_to_mean,
            )
        finally:
            self._release_file(lock)

        return output_path

    def _append_summary(self, df: pd.DataFrame):
        sum_lock = self._lock_file(self.summary_path)
        try:
            csv_path = self.summary_path
            header = not os.path.isfile(csv_path)
            df.to_csv(csv_path, mode="a", index=False, header=header)

            root, _ = os.path.splitext(csv_path)
            parquet_path = f"{root}.parquet"
            pq.write_table(
                pa.Table.from_pandas(pd.read_csv(csv_path), preserve_index=False),
                parquet_path,
                compression="snappy",
            )
        finally:
            self._release_file(sum_lock)

    def _lock_file(self, path: str) -> FileLock:
        lock = FileLock(path + ".lock")
        lock.acquire()
        return lock

    def _release_file(self, lock: FileLock):
        lock.release()
        if os.path.exists(lock.lock_file):
            os.remove(lock.lock_file)