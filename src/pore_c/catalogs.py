from intake.catalog.local import YAMLFileCatalog
from pathlib import Path
import yaml


class basePoreCCatalog(YAMLFileCatalog):
    _suffix_map ={
    }
    @classmethod
    def generate_paths(cls, prefix):
        res = {}
        exists = []
        for key, val in cls._suffix_map.items():
            res[key] = Path(prefix + val)
            if res[key].exists():
                exists.append(key)
        return res, exists


class VirtualDigestCatalog(basePoreCCatalog):
    name = 'pore_c_virtual_digest'
    description = "An intake catalog file for a virtual digest of a reference genome"

    _suffix_map ={
            'catalog': '.catalog.yaml',
            'fragments': ".fragments.parquet",
            'digest_stats': ".digest_stats.csv",
    }

    @classmethod
    def create(cls, file_paths, refgenome_catalog, digest_type, digest_param, num_fragments):
        catalog_path = file_paths.pop('catalog')
        catalog_data = {
            "name": cls.name,
            "description": cls.description,
            "driver": 'pore_c.catalogs.VirtualDigestCatalog',
            "metadata": {
                "digest_type": digest_type,
                "digest_param": digest_param,
                "num_fragments": num_fragments,
            },
            "sources": {
                "refgenome": {
                    "args": {"path": str(refgenome_catalog.resolve())},
                    "driver": "intake.catalog.local.YAMLFileCatalog",
                    "description": ReferenceGenomeCatalog.description
                }
            },
        }

        for key, val in file_paths.items():
            driver = val.suffix.replace('.', '')
            assert(driver in ['parquet', 'csv'])
            catalog_data['sources'][key] = {
                'driver': driver,
                'args': {'urlpath': str(val.resolve())}
            }
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat

    def __str__(self):
        return "<VirtualDigestCatalog digest_type={} digest_param:{} num_fragments:{} path:{}>".format(
            self.metadata['digest_type'], self.metadata['digest_param'], self.metadata['num_fragments'], self.path
        )






class ReferenceGenomeCatalog(basePoreCCatalog):
    name = 'pore_c_reference_genome'
    description = "An intake catalog file for a reference genome"

    _suffix_map ={
            'catalog': '.catalog.yaml',
            'chrom_metadata': ".metadata.csv"
    }

    @classmethod
    def create(cls, catalog_path, fasta_path, metadata_csv, chrom_lengths, genome_id):
        catalog_data = {
            "name": cls.name,
            "description": cls.description,
            "driver": 'pore_c.catalogs.ReferenceGenomeCatalog',
            "metadata": {
                "chrom_lengths": chrom_lengths,
                "genome_id": genome_id
            },
            "sources": {
                "fasta": {"driver": "pore_c.datasources.IndexedFasta", "args": {"urlpath": "{}".format(fasta_path.resolve())}},
                "chrom_metadata": {"driver": "csv", "args": {"urlpath": "{{ CATALOG_DIR }}/" + str(metadata_csv.name)}},
            },
        }
        with catalog_path.open("w") as fh:
            fh.write(yaml.dump(catalog_data, default_flow_style=False, sort_keys=False))
        cat = cls(str(catalog_path))
        return cat

    @property
    def chrom_lengths(self):
        return self.metadata['chrom_lengths']

    @property
    def chrom_order(self):
        return list(self.metadata['chrom_lengths'].keys())

    @property
    def chrom_dtype(self):
        return pd.CategoricalDtype(self.chrom_order, ordered=True)

    @property
    def genome_size(self):
        return sum(self.metadata['chrom_lengths'].values())

    @property
    def genome_id(self):
        return self.metadata['genome_id']

    def __str__(self):
        return "<ReferenceGenomeCatalog genome_id={} genome_size={:,} num_chroms={} chroms={}..{}>".format(
            self.genome_id, self.genome_size, len(self.chrom_lengths), ",".join(self.chrom_order[:3]), self.chrom_order[-1])




