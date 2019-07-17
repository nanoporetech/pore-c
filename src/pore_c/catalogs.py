from pathlib import Path

import yaml
import pandas as pd
from intake.catalog.local import YAMLFileCatalog


class basePoreCCatalog(YAMLFileCatalog):
    _suffix_map = {}

    @classmethod
    def generate_paths(cls, prefix):
        res = {}
        exists = []
        for key, val in cls._suffix_map.items():
            res[key] = Path(prefix + val)
            if res[key].exists():
                exists.append(key)
        return res, exists


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
            catalog_data["sources"][key] = {"driver": driver, "args": {"urlpath": str(val.resolve())}}
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
            catalog_data["sources"][key] = {"driver": driver, "args": {"urlpath": str(val.resolve())}}
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

    _suffix_map = {"catalog": ".catalog.yaml", "fragments": ".fragments.parquet", "digest_stats": ".digest_stats.csv"}

    @classmethod
    def create(cls, file_paths, refgenome_catalog, digest_type, digest_param, num_fragments):
        catalog_path = file_paths.pop("catalog")
        catalog_data = {
            "name": cls.name,
            "description": cls.description,
            "driver": "pore_c.catalogs.VirtualDigestCatalog",
            "metadata": {"digest_type": digest_type, "digest_param": digest_param, "num_fragments": num_fragments},
            "sources": {
                "refgenome": {
                    "args": {"path": str(refgenome_catalog.resolve())},
                    "driver": "intake.catalog.local.YAMLFileCatalog",
                    "description": ReferenceGenomeCatalog.description,
                }
            },
        }

        for key, val in file_paths.items():
            driver = val.suffix.replace(".", "")
            assert driver in ["parquet", "csv"]
            catalog_data["sources"][key] = {"driver": driver, "args": {"urlpath": str(val.resolve())}}
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat

    def __str__(self):
        return "<VirtualDigestCatalog digest_type={} digest_param:{} num_fragments:{} path:{}>".format(
            self.metadata["digest_type"], self.metadata["digest_param"], self.metadata["num_fragments"], self.path
        )


class ReferenceGenomeCatalog(basePoreCCatalog):
    name = "pore_c_reference_genome"
    description = "An intake catalog file for a reference genome"

    _suffix_map = {"catalog": ".catalog.yaml", "chromsizes": ".chromsizes", "chrom_metadata": ".metadata.csv"}

    @classmethod
    def create(cls, catalog_path, fasta_path, metadata_csv, chrom_lengths, chromsizes, genome_id):
        catalog_data = {
            "name": cls.name,
            "description": cls.description,
            "driver": "pore_c.catalogs.ReferenceGenomeCatalog",
            "metadata": {"chrom_lengths": chrom_lengths, "genome_id": genome_id},
            "sources": {
                "fasta": {
                    "driver": "pore_c.datasources.IndexedFasta",
                    "args": {"urlpath": "{}".format(fasta_path.resolve())},
                },
                "chrom_metadata": {"driver": "csv", "args": {"urlpath": "{{ CATALOG_DIR }}/" + str(metadata_csv.name)}},
                "chromsizes": {"driver": "csv", "args": {"urlpath": "{{ CATALOG_DIR }}/" + str(chromsizes.name)}},
            },
        }
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
