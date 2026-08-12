# H3SummaryStager.py (Python 3.9 compatible)

import logging
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Union, List

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from filelock import FileLock

from . import TilePathManager
from .H3GridSummaryGenerator import H3GridSummaryGenerator


PathLike = Union[str, Path]


class H3SummaryStager:
    """
    Runs H3GridSummaryGenerator over all input vectors and writes one
    H3 summary GeoPackage per input file under a configured base_dir.
    Optionally records a CSV/Parquet run summary.
    """

    def __init__(
        self,
        config,
        tiles: Optional[TilePathManager] = None,
        out_base_dir: str = "h3",
        out_ext: str = ".gpkg",
        summary_filename: str = "h3_summary.csv",
        generator: Optional[H3GridSummaryGenerator] = None,
    ):
        self.logger = logging.getLogger(__name__)
        self.config = config
        self.tiles = tiles
        self.out_base_dir = out_base_dir
        self.out_ext = out_ext

        if self.tiles is None:
            raise ValueError("tiles (TilePathManager) is required for workflow integration")

        h3_root = Path(self.tiles.base_dirs["h3"]["path"])
        h3_root.mkdir(parents=True, exist_ok=True)
        self.summary_path = str(h3_root / summary_filename)

        self.gen = generator or H3GridSummaryGenerator(
            tiles=self.tiles,
            out_base_dir=out_base_dir,
            logger=self.logger,
        )

    def _output_path_for_input(self, input_path: PathLike, h3_res: int) -> Path:
        input_path = Path(input_path)

        if self.config.is_stager_enabled():
            input_root = Path(self.tiles.base_dirs["staged"]["path"])
        else:
            input_root = Path(self.tiles.base_dirs["input"]["path"])

        try:
            rel = input_path.relative_to(input_root)
        except Exception:
            rel = Path(input_path.name)

        stem = rel.stem
        parent = rel.parent
        out_root = Path(self.tiles.base_dirs["h3"]["path"])
        out_dir = out_root / parent
        out_dir.mkdir(parents=True, exist_ok=True)

        return out_dir / f"{stem}_h3r{h3_res}{self.out_ext}"

    def stage_all(
        self,
        h3_res: List[int],
        attr_to_sum: Optional[List[str]] = None,
        attr_to_mean: Optional[List[str]] = None,
        land_polygons_path: Optional[PathLike] = None,
        area_epsg: Optional[int] = None,
    ) -> None:
        overall_start = datetime.now()

        if self.config.is_stager_enabled():
            input_paths = self.tiles.get_filenames_from_dir("staged")
        else:
            input_paths = self.tiles.get_filenames_from_dir("input")
        n = len(input_paths)

        if n == 0:
            self.logger.error("No vector files found for H3 staging.")
            return

        out_root = Path(self.tiles.base_dirs["h3"]["path"])
        out_root.mkdir(parents=True, exist_ok=True)
        
        # build a dictionary of paths for every requested resolution
        master_outputs = {}
        for res in h3_res:
            out_path = out_root / f"summary_h3r{res}{self.out_ext}"
            if out_path.exists():
                self.logger.info("Removing existing H3 file to start fresh: %s", out_path)
                out_path.unlink()
            master_outputs[res] = out_path

        self.logger.info("Begin H3 staging %s input vector files.", n)

        rows = []
        for p in input_paths:
            start = datetime.now()
            out_dict = None
            ok = False
            err = None
            try:
                out_dict = self.stage(
                    path=p,
                    h3_res=h3_res,
                    attr_to_sum=attr_to_sum,
                    attr_to_mean=attr_to_mean,
                    land_polygons_path=land_polygons_path,
                    area_epsg=area_epsg,
                    output_paths=master_outputs
                )
                ok = True
            except Exception as e:
                err = repr(e)
                self.logger.exception("Failed to stage H3 summary for %s", p)

            rows.append(
                {
                    "input_path": str(p),
                    "output_paths": str(out_dict) if out_dict else None, 
                    "h3_res": str(h3_res),
                    "ok": ok,
                    "seconds": (datetime.now() - start).total_seconds(),
                    "error": err,
                    "datetime": datetime.now().isoformat(),
                }
            )

        df = pd.DataFrame(rows)
        self._append_summary(df)

        total = datetime.now() - overall_start
        self.logger.info("H3-staged %s files in %s (%s per file).", n, total, total / max(n, 1))

    def stage(
        self,
        path: PathLike,
        h3_res: List[int],
        attr_to_sum: Optional[List[str]] = None,
        attr_to_mean: Optional[List[str]] = None,
        land_polygons_path: Optional[PathLike] = None,
        area_epsg: Optional[int] = None,
        output_paths: Optional[dict] = None,
    ) -> dict:
        
        if output_paths is None:
            out_root = Path(self.tiles.base_dirs["h3"]["path"])
            out_root.mkdir(parents=True, exist_ok=True)
            output_paths = {
                res: out_root / f"summary_h3r{res}{self.out_ext}" 
                for res in h3_res
            }

        locks = []
        # Acquire locks for all resolution files
        for out_path in output_paths.values():
            locks.append(self._lock_file(str(out_path)))
            
        try:
            self.gen.build_h3_summary(
                input_path=path,
                output_paths=output_paths,
                h3_res=h3_res,
                land_polygons_path=land_polygons_path,
                area_epsg=area_epsg,
                attr_to_sum=attr_to_sum,
                attr_to_mean=attr_to_mean,
            )
        finally:
            # Ensure all locks are released even if a failure occurs
            for lock in locks:
                self._release_file(lock)
                
        return output_paths

    def _append_summary(self, df: pd.DataFrame) -> None:
        # lock summary so concurrent stage_all() runs don't corrupt the log
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

    def _release_file(self, lock: FileLock) -> None:
        lock.release()
        if os.path.exists(lock.lock_file):
            os.remove(lock.lock_file)
