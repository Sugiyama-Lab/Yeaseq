from json import load

from .yeaseq_config import WorkDir


class OrthoGenes(object):
    FileName_OrthoGeneTable = 'OrthoGeneTable-Rhind_Science_11.json'

    def __init__(self):
        self.ortho_table_path, self.ortho_table = self._load_ortho_gene_table()

    def _load_ortho_gene_table(self):
        ortho_table_path = WorkDir.get_external_data_path(self.FileName_OrthoGeneTable)
        with open(ortho_table_path, 'r') as f:
            ortho_table = load(f)
        return ortho_table_path, ortho_table

    def find_gene(self, gene: str):
        record_idx = self.ortho_table.get(gene)
        if record_idx is None:
            record_idx = self.ortho_table.get(gene.lower())
        if record_idx is None:
            record_idx = self.ortho_table.get(gene.upper())
        if record_idx is None:
            return None

        record = self.ortho_table['Data'][record_idx]
        return record

    @staticmethod
    def fill_shown_text(gene, record):
        texts = [
            '\nThe orthologous information of {}:\n\n'.format(gene),

            ('  S. pombe: {} with gene name {}\n'.format(record[1], record[0])
             if record[0] != ''
             else '  S. pombe: {}\n'.format(record[1])),

            '  S. japonicus: {}\n'.format(record[2]),
            '  S. cryophilus: {}\n'.format(record[3]),
            '  S. octosporus: {}\n'.format(record[4]),

            '\nThe gene description of S. pombe: {}.'.format(record[5])

        ]
        return texts
