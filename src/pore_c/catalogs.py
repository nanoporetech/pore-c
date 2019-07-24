from logging import getLogger
from pathlib import Path

import pandas as pd
import yaml
from intake.catalog.local import YAMLFileCatalog
from intake.source.base import DataSource

__all__ = [
    "RawReadCatalog",
    "PairsFileCatalog",
    "AlignmentDfCatalog",
    "VirtualDigestCatalog",
    "ReferenceGenomeCatalog",
]

logger = getLogger(__name__)


class basePoreCCatalog(YAMLFileCatalog):
    _suffix_map = {}

    @classmethod
    def generate_paths(cls, prefix):
        res = {}
        exists = []
        for key, val in cls._suffix_map.items():
            res[key] = Path(prefix + val)
            logger.debug(
                "Data for {} {} will be saved to: {}".format(cls.__name__, key, val)
            )
            if res[key].exists():
                exists.append(key)
        if exists:
            for file_id in exists:
                logger.error(
                    "Output file already exists for {}: {}".format(
                        file_id, res[file_id]
                    )
                )
            raise IOError("Please remove output files before continuing")
        return res

    @property
    def user_metadata(self):
        return self.metadata.get("user_metadata", {})

    @property
    def porec_metadata(self):
        md = self.metadata.copy()
        if "user_metadata" in md:
            md.pop("user_metadata")
        return md

    @staticmethod
    def create_catalog_dict(klass, file_paths, metadata, user_metadata, drivers=None):
        if drivers is None:
            drivers = {}
        if not metadata:
            metadata = {}
        if user_metadata:
            metadata["user_metadata"] = user_metadata
        catalog_data = {
            "name": klass.name,
            "description": klass.description,
            "driver": "pore_c.catalogs.{}".format(klass.__name__),
            "metadata": metadata,
            "sources": {},
        }
        driver_lookup = {
            ".csv": "csv",
            ".parquet": "parquet",
            ".fq.gz": "pore_c.datasources.Fastq",
            ".fasta.gz": "pore_c.datasources.IndexedFasta",
            ".fq": "pore_c.datasources.Fastq",
            ".pairs.gz": "pore_c.datasources.IndexedPairFile",
            ".catalog.yaml": "intake.catalog.local.YAMLFileCatalog"
        }
        catalog_path = file_paths["catalog"]
        for key, val in file_paths.items():
            if key == "catalog":
                pass
            elif isinstance(val, DataSource):
                raise NotImplementedError
                catalog_data['sources'][key] = val._yaml()
            elif isinstance(val, Path):
                driver = drivers.get(key, None)
                if driver is None:
                    for suffix, _driver in driver_lookup.items():
                        if val.name.endswith(suffix):
                            driver = _driver
                            break
                if driver is None:
                    logger.warning(f"No driver found found for {key}: {val}, skipping.")
                    continue
                try:
                    urlpath = "{{ CATALOG_DIR }}/" + str(
                        val.relative_to(catalog_path.parent)
                    )
                except ValueError:
                    urlpath = str(val.resolve())
                catalog_data["sources"][key] = {
                    "driver": driver,
                    "args": {"urlpath": urlpath},
                }
            else:
                raise ValueError(val)
        return catalog_data


class RawReadCatalog(basePoreCCatalog):
    name = "pore_c_raw_reads"
    description = "A catalog of data associated with raw reads"

    _suffix_map = {
        "catalog": ".catalog.yaml",
        "pass_fastq": ".pass.fq.gz",
        "fail_fastq": ".fail.fq.gz",
        "read_metadata": ".read_metadata.parquet",
        "summary": ".summary.csv",
    }

    @classmethod
    def create(cls, file_paths, *args, **kwds):
        catalog_data = basePoreCCatalog.create_catalog_dict(
            cls, file_paths, *args, **kwds
        )
        catalog_path = file_paths["catalog"]
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat

    def __str__(self):
        return "<RawReadCatalog reads={num_sequences} bases={total_bases:,}>".format(
            **self.metadata["pass"]
        )


class PairsFileCatalog(basePoreCCatalog):
    name = "pore_c_pairs"
    description = "An intake catalog file for a pairs format file"

    _suffix_map = {"catalog": ".catalog.yaml", "pairs": ".pairs.gz"}

    @classmethod
    def create(cls, file_paths, alignment_df_catalog):
        catalog_path = file_paths.pop("catalog")
        catalog_data = {
            "name": cls.name,
            "description": cls.description,
            "sources": {
                "source_alignments": {
                    "args": {"path": str(alignment_df_catalog.resolve())},
                    "driver": "intake.catalog.local.YAMLFileCatalog",
                    "description": AlignmentDfCatalog.description,
                }
            },
        }
        for key, val in file_paths.items():
            if key == "pairs":
                driver = "pore_c.datasources.IndexedPairFile"
            else:
                driver = val.suffix.replace(".", "")
                assert driver in ["parquet", "csv"]
            catalog_data["sources"][key] = {
                "driver": driver,
                "args": {"urlpath": str(val.resolve())},
            }
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat


class AlignmentDfCatalog(basePoreCCatalog):
    name = "pore_c_alignment_df"
    description = "An intake catalog file for the annotated alignments dataframe"

    _suffix_map = {
        "catalog": ".catalog.yaml",
        "alignment": ".alignment.parquet",
        "read": ".read.parquet",
        "overlap": ".overlap.parquet",
        "alignment_summary": ".alignment_summary.csv",
        "read_summary": ".read_summary.csv",
    }

    @classmethod
    def create(cls, file_paths, input_bam, virtual_digest_catalog, final_stats):
        catalog_path = file_paths.pop("catalog")
        catalog_data = {
            "name": "pore_c_parsed_alignment_files",
            "description": "Output files of pore-c tools alignment filtering",
            "sources": {
                "virtual_digest": {
                    "args": {"path": str(virtual_digest_catalog.resolve())},
                    "driver": "intake.catalog.local.YAMLFileCatalog",
                    "description": VirtualDigestCatalog.description,
                },
                "source_bam": {
                    "args": {"urlpath": str(input_bam.resolve())},
                    "driver": "pore_c.datasources.NameSortedBamSource",
                    "description": "The input bam file",
                },
            },
            "metadata": final_stats,
        }
        for key, val in file_paths.items():
            driver = val.suffix.replace(".", "")
            assert driver in ["parquet", "csv"]
            catalog_data["sources"][key] = {
                "driver": driver,
                "args": {"urlpath": str(val.resolve())},
            }
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat

    def __str__(self):
        return "<AlignmentDfCatalog>"
        # ={} digest_param:{} num_fragments:{} path:{}>".format(
        #    self.metadata['digest_type'], self.metadata['digest_param'], self.metadata['num_fragments'], self.path
        # )


class VirtualDigestCatalog(basePoreCCatalog):
    name = "pore_c_virtual_digest"
    description = "An intake catalog file for a virtual digest of a reference genome"

    _suffix_map = {
        "catalog": ".catalog.yaml",
        "fragments": ".fragments.parquet",
        "digest_stats": ".digest_stats.csv",
    }

    @classmethod
    def create(cls, file_paths, *args, **kwds):
        catalog_data = basePoreCCatalog.create_catalog_dict(
            cls, file_paths, *args, **kwds
        )
        catalog_path = file_paths["catalog"]
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))

        #catalog_path = file_paths.pop("catalog")
        #catalog_data = {
        #    "name": cls.name,
        #    "description": cls.description,
        #    "driver": "pore_c.catalogs.VirtualDigestCatalog",
        #    "metadata": {
        #        "digest_type": digest_type,
        #        "digest_param": digest_param,
        #        "num_fragments": num_fragments,
        #    },
        #    "sources": {
        #        "refgenome": {
        #            "args": {"path": str(refgenome_catalog.resolve())},
        #            "driver": "intake.catalog.local.YAMLFileCatalog",
        #            "description": ReferenceGenomeCatalog.description,
        #        }
        #    },
        #}

        #for key, val in file_paths.items():
        #    driver = val.suffix.replace(".", "")
        #    assert driver in ["parquet", "csv"]
        #    catalog_data["sources"][key] = {
        #        "driver": driver,
        #        "args": {"urlpath": str(val.resolve())},
        #    }
        #with catalog_path.open("w") as fh:
        #    fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        #cat = cls(str(catalog_path))
        #return cat

    def __str__(self):
        return "<VirtualDigestCatalog digest_type={} digest_param:{} num_fragments:{} path:{}>".format(
            self.metadata["digest_type"],
            self.metadata["digest_param"],
            self.metadata["num_fragments"],
            self.path,
        )


class ReferenceGenomeCatalog(basePoreCCatalog):
    name = "pore_c_reference_genome"
    description = "An intake catalog file for a reference genome"

    _suffix_map = {
        "catalog": ".catalog.yaml",
        "chromsizes": ".chromsizes",
        "chrom_metadata": ".metadata.csv",
    }

    @classmethod
    def create(cls, file_paths, *args, **kwds):
        catalog_data = basePoreCCatalog.create_catalog_dict(
                cls, file_paths, *args, **kwds, drivers={'chromsizes': 'csv'}
        )
        catalog_path = file_paths["catalog"]
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat

    @property
    def chrom_lengths(self):
        return self.metadata["chrom_lengths"]

    @property
    def chrom_order(self):
        return list(self.metadata["chrom_lengths"].keys())

    @property
    def chrom_dtype(self):
        return pd.CategoricalDtype(self.chrom_order, ordered=True)

    @property
    def genome_size(self):
        return sum(self.metadata["chrom_lengths"].values())

    @property
    def genome_id(self):
        return self.metadata["genome_id"]

    def __str__(self):
        return "<ReferenceGenomeCatalog genome_id={} genome_size={:,} num_chroms={} chroms={}..{}>".format(
            self.genome_id,
            self.genome_size,
            len(self.chrom_lengths),
            ",".join(self.chrom_order[:3]),
            self.chrom_order[-1],
        )
