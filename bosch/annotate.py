import csv
import os
from typing import Generator, List, Tuple

from tqdm import tqdm

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pyhmmer.easel import Alphabet, TextSequence, DigitalSequence
from pyhmmer.plan7 import HMMFile, OptimizedProfile
from pyhmmer.hmmer import hmmscan


def _get_protein_sequences(record: SeqRecord) -> List[DigitalSequence]:
    '''...

    Args:
        record (SeqRecord): ...

    Returns:
        List[DigitalSequence]: ...
    '''
    alphabet = Alphabet.amino()
    protein_sequences = []

    for feature in record.features:
        if feature.type == 'CDS':
            if 'translation' in feature.qualifiers:
                if 'locus_tag' in feature.qualifiers:
                    seq_name = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    seq_name = feature.qualifiers['gene'][0]
                else:
                    continue

                seq = feature.qualifiers['translation'][0]
                seq = TextSequence(name=seq_name.encode(), sequence=seq)
                
                protein_sequences.append(seq.digitize(alphabet=alphabet))
    
    return protein_sequences


def _assign_pfam_family(queries: List[DigitalSequence],
                        profiles: List[OptimizedProfile],
                        cpus: int,
                        E: float,
                        verbose: bool) -> Generator[Tuple[str, str, str], None, None]:
    '''...
    
    Args:
        queries (List[DigitalSequence]): ...
        profiles (List[OptimizedProfile]): ...
        cpus (int): ...
        E (float): ...
    Yields:
        Tuple[str, str, str]: ...
    '''
    for hits in tqdm(hmmscan(queries, profiles, cpus=cpus, E=E),
                     total=len(queries), disable=not verbose):
        query_name = hits.query_name.decode()

        if hits:
            pfam_acc = hits[0].accession.decode()
            pfam_desc = hits[0].description.decode()

            yield query_name, pfam_acc, pfam_desc
        else:
            yield query_name, 'PFAM_UNK', ''


def run(src: os.PathLike, dst: os.PathLike, cpus: int, E: float, verbose: bool = True) -> None:
    '''...

    Args:
        src (os.PathLike): ...
        dst (os.PathLike): ...
        cpus (int): ...
        E (float): ...
    '''
    PFAM = os.path.join(os.path.dirname(__file__), 'Pfam', 'Pfam-A.hmm')

    if os.path.exists(PFAM):
        with HMMFile(PFAM) as file:
            profiles = list(file.optimized_profiles())
    else:
        raise FileNotFoundError

    for record in SeqIO.parse(src, 'genbank'):
        queries = _get_protein_sequences(record)

        with open(os.path.join(dst, record.id + '.tsv'), mode='w') as file:
            writer = csv.writer(file, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(('locus_tag', 'pfam_acc', 'pfam_desc'))

            for hit in _assign_pfam_family(queries, profiles, cpus, E, verbose=verbose):
                writer.writerow(hit)
