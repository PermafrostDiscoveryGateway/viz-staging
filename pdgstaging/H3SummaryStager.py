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
    """Executes H3 grid summarization across input vector datasets.

    This class runs the H3GridSummaryGenerator over all provided input 
    vectors and writes one H3 summary GeoPackage per input file into a 
    configured base directory. It can optionally record a run summary in 
    CSV or Parquet format for downstream tracking.
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
        """Initializes the H3SummaryStager.

        Args:
            config: The configuration manager or dictionary containing global settings.
            tiles (Optional[TilePathManager], optional): Manager for resolving 
                input and output tile paths. Defaults to None.
            out_base_dir (str, optional): The base directory for writing output 
                GeoPackage files. Defaults to "h3".
            out_ext (str, optional): The file extension for the generated 
                summaries. Defaults to ".gpkg".
            summary_filename (str, optional): The name of the run summary file 
                used for tracking. Defaults to "h3_summary.csv".
            generator (Optional[H3GridSummaryGenerator], optional): The generator 
                instance used to process the grid summaries. Defaults to None.
        """
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
        """Determines the target output file path for a given input path and H3 resolution.

        Args:
            input_path (PathLike): The file path to the input dataset.
            h3_res (int): The H3 resolution level used for file naming.

        Returns:
            Path: The fully resolved, absolute destination path for the output file.
        """
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
        """Executes full H3 staging and aggregation across all input vector files.

        Iterates over input vector files (staged or raw, depending on configuration),
        generates intermediate H3 summary chunks for each requested resolution, and
        combines those chunks into unified output GeoPackage files. Execution metrics
        and errors for each file are tracked and logged to a run summary file.

        Args:
            h3_res (List[int]): A list of H3 resolution levels to generate summaries for.
            attr_to_sum (Optional[List[str]], optional): List of column names to aggregate 
                using summation. Defaults to None.
            attr_to_mean (Optional[List[str]], optional): List of column names to aggregate 
                using mean calculations. Defaults to None.
            land_polygons_path (Optional[PathLike], optional): File path to a vector 
                dataset of land boundaries used for land area calculations. Defaults to None.
            area_epsg (Optional[int], optional): EPSG code for the equal-area CRS used 
                in area calculations. Defaults to None.
        """
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
        
        self.logger.info("All chunks generated. Combining into final GeoPackages...")
        self.gen.combine_h3_summaries(
            output_paths=master_outputs,
            h3_res=h3_res,
            land_polygons_path=land_polygons_path,
            area_epsg=area_epsg,
            attr_to_sum=attr_to_sum,
            attr_to_mean=attr_to_mean
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
        """Processes a single input vector file and builds its H3 summary chunks.

        Resolves the input path, ensures the intermediate output structure exists,
        and delegates execution to the generator to create intermediate Parquet
        chunks for each requested H3 resolution. All resolutions are done at the same
        time to ensure each input file is only opened once.

        Args:
            path (PathLike): File path to the input vector dataset to process.
            h3_res (List[int]): A list of H3 resolution levels to generate 
                summaries for.
            attr_to_sum (Optional[List[str]], optional): List of column names to 
                aggregate using summation. Defaults to None.
            attr_to_mean (Optional[List[str]], optional): List of column names to 
                aggregate using mean calculations. Defaults to None.
            land_polygons_path (Optional[PathLike], optional): Path to the vector 
                dataset of land boundaries. Defaults to None.
            area_epsg (Optional[int], optional): EPSG code for the equal-area CRS 
                used in area calculations. Defaults to None.
            output_paths (Optional[dict], optional): A dictionary mapping H3 
                resolutions (int) to their target final output paths (PathLike). 
                Defaults to None.

        Returns:
            dict: A dictionary mapping each H3 resolution to its corresponding 
                generated output path.
        """
        
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
        """Appends execution metrics to a CSV log and updates a Parquet summary file.

        Args:
            df (pd.DataFrame): DataFrame containing summary rows and execution 
                metrics to log.
        """
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
        """Acquires a file lock to prevent concurrent process race conditions.

        Args:
            path (str): The file path for which to create the lock.

        Returns:
            FileLock: An acquired `FileLock` object managing the lock lifecycle.
        """
        lock = FileLock(path + ".lock")
        lock.acquire()
        return lock

    def _release_file(self, lock: FileLock) -> None:
        """Releases an acquired file lock and removes the associated lock file.

        Args:
            lock (FileLock): The `FileLock` instance to release and clean up.
        """
        lock.release()
        if os.path.exists(lock.lock_file):
            os.remove(lock.lock_file)
