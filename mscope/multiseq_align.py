import os
import binascii
import numpy
import Levenshtein
import sys

from Bio import SeqIO,Align
from mscope.util import MotifPlacement
import pandas as pd
import pyabpoa

def generate_hex_chr_list():
    usable_hex = list(range(1,256))
    usable_hex.remove(60)
    usable_hex.remove(61)
    usable_hex.remove(62)
    usable_hex.remove(45)
    usable_hex.remove(32)
    usable_hex.remove(10)
    usable_hex.remove(13)
    usable_hex = ['{:02x}'.format(x) for x in usable_hex]
    return usable_hex

def assign_chr_to_motif(all_unique_motif, printable_asii):
    if len(all_unique_motif) > len(printable_asii):
        return None, "too many motifs"        
    else:
        printable_asii_to_include = printable_asii[:len(all_unique_motif)]
        motif_dict = dict(zip(all_unique_motif, printable_asii_to_include))
        return motif_dict, printable_asii_to_include

def hex_to_bytes(hex_string):
    return binascii.unhexlify(hex_string)

def normalize_edit_distance(distance, seq_1_length, seq_2_length):
    if max(seq_1_length, seq_2_length) != 0:
        return distance / max(seq_1_length, seq_2_length)
    else:
        return 0


class Aligner:
    """Base class for aligners.
    """
    def __init__(self, cfg):
        self.filename_prefix = cfg.output
        self.cfg = cfg
        self.printable_asii_list = generate_hex_chr_list()

    def _calc_sub_matrix(self, hex_used, motif_dict):
        """Calculate the substitution score matrix for the motifs
           
        :param hex_used: list of characters used in the alignment
        :param motif_dict: dict, mapping of motifs to characters

        :returns: numpy.array, the substitution score matrix for the motifs
        """
        sub_matrix = numpy.zeros((len(hex_used), len(hex_used)))
        rev_motif_dict = {y: x for x, y in motif_dict.items()}
        paligner = Align.PairwiseAligner()
        paligner.mode = 'global'
        paligner.match_score = 1
        paligner.mismatch_score = -1
        paligner.open_gap_score = -2
        paligner.extend_gap_score = -1

        for i in range(len(hex_used)):
            for j in range(i, len(hex_used)):
                seq1 = rev_motif_dict[hex_used[i]]
                seq2 = rev_motif_dict[hex_used[j]]
                sub_matrix[i,j] = paligner.score(seq1, seq2)
                sub_matrix[j,i] = sub_matrix[i,j]

        return sub_matrix

    def _calculate_edit_distance(self, result, seq_ids):
        """Calculate the edit distance between the sequences

        :param result: numpy.array, the multiple sequence alignment result. nseq * sequence_length. Each cell contains a motif or a base
        :param seq_ids: list of sequence ids


        :returns: pandas.DataFrame, the edit distance between the sequences
        """
        nseq, sequence_length = result.shape

        result_filled = result != ''

        dist = numpy.zeros((nseq, nseq))
        for i in range(len(seq_ids)):
            for j in range(i+1, len(seq_ids)):
                fil = result_filled[i,:] | result_filled[j,:]
                common = fil.sum()
                tdist = (result[i,:] !=  result[j,:]).sum()
                dist[i, j] = normalize_edit_distance(tdist, common, common)
                dist[j,i] = dist[i,j]
        return pd.DataFrame(dist, index=seq_ids, columns=seq_ids)                
       

    def _result_to_grouped_positions(self, result, seq_ids):
        """Convert the multiple sequence motif alignment result to per sequence a list of motifplacement objects.

        :param result: numpy.array, the multiple sequence alignment result. nseq * sequence_length. Each cell contains a motif or a base
        """
        nseq, sequence_length = result.shape
        unique_motifs = [set(result[:,i]) for i in range(sequence_length)]
        max_length = [0] + [max([len(e) for e in unique_motifs[i]]) for i in range(sequence_length)]
        pos = numpy.cumsum(max_length)
        start = pos[:-1]
        stop = pos[1:] 

        grouped_result = {}
        for i in range(len(seq_ids)):
            motifs = []
            for j in range(sequence_length):
                seq = result[i,j]
                if not seq:
                    continue

                if not (motifs and motifs[-1].attempt_add(seq, start=start[j])):
                    motifs.append(MotifPlacement(seq, start[j]))

            key = seq_ids[i]
            grouped_result[key] = motifs

        return grouped_result, stop[-1]

class LevenshteinAligner(Aligner):
    """Aligner that uses the Levenshtein distance to calculate the edit distance between the sequences.
       No multiple sequence alignment is performed.
    """
    def __init__(self, cfg, all_seq_dict):
        Aligner.__init__(self, cfg)
        self.all_seq_dict = all_seq_dict


    def run(self, grouped_positions, all_unique_motifs):
        sys.stderr.write("No multiple sequence alignment, using Levenshtein edit distance ...\n"); sys.stderr.flush()
        motif_dict, hex_used = assign_chr_to_motif(all_unique_motifs, self.printable_asii_list)
        
        sequence_lengths = {seq: len(self.all_seq_dict[seq]) for seq in self.all_seq_dict}

        ngrouped_positions = {}
        for sequence, motifs in grouped_positions.items():
            
            seq = ''.join([''.join([motif_dict[mp.motif]] * mp.count) if not mp.singlebase else ''.join([motif_dict[mc] for mc in mp.motif])   for mp in motifs])
            if seq:
                ngrouped_positions[sequence] = hex_to_bytes(seq)
            else:
                ngrouped_positions[sequence] = b" "
        
        sequences = list(ngrouped_positions.keys())
        dist = numpy.zeros((len(sequences), len(sequences)))
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                seq1 = ngrouped_positions[sequences[i]]
                seq2 = ngrouped_positions[sequences[j]]
                tdist = Levenshtein.distance(seq1, seq2)
                dist[i, j] = normalize_edit_distance(tdist, len(seq1), len(seq2))
                dist[j,i] = dist[i,j]

        return grouped_positions, sequence_lengths, pd.DataFrame(dist, index=sequences, columns=sequences) 

class MAFFTAligner(Aligner):
    """Aligner that uses MAFFT to perform the multiple sequence alignment and edit distance calculation.
    """
    def __init__(self, cfg):
        Aligner.__init__(self, cfg)
        self.hex_fasta_file = self.filename_prefix + "_motif_in_hex" + ".hex"
        self.msa_sub_matrix = self.filename_prefix + "_msa_sub_matrix" + ".txt"

    def run(self, grouped_positions, all_unique_motifs):
        sys.stderr.write("Aligning sequences using MAFFT aligner...\n"); sys.stderr.flush()
        motif_dict, hex_used = assign_chr_to_motif(all_unique_motifs, self.printable_asii_list)

        sub_matrix = self._calc_sub_matrix(hex_used, motif_dict)
        self.write_sub_matrix(sub_matrix, hex_used, motif_dict)
        
        self.write_seq_in_hex_chr(grouped_positions, motif_dict)
        msa_result, msa_keys = self._msa_with_characters(motif_dict)
        
        edit_distance = self._calculate_edit_distance(msa_result, msa_keys)

        ngrouped_positions, align_length = self._result_to_grouped_positions(msa_result, msa_keys)

        sequence_lengths = {key: align_length for key, value in grouped_positions.items()}
        return ngrouped_positions, sequence_lengths, edit_distance


    def write_seq_in_hex_chr(self, grouped_positions, motif_dict):
        with open (self.hex_fasta_file, 'w') as hex_fasta:
            for sequence, motifs in grouped_positions.items():
                seq = ' '.join([' '.join([motif_dict[mp.motif]] * mp.count) if not mp.singlebase else ' '.join([motif_dict[mc] for mc in mp.motif])   for mp in motifs])
                hex_fasta.write(">" + sequence + "\n")
                hex_fasta.write(seq + "\n")
        

    def write_sub_matrix(self, sub_matrix, hex_used, motif_dict):
        rev_motif_dict = {y: x for x, y in motif_dict.items()}
        with open(self.msa_sub_matrix, 'w') as sub_matrix_file:
            for i in range(len(hex_used)):
                for j in range(i, len(hex_used)):
                    #hex_used has integer values
                    line = f"0x{hex_used[i]} 0x{hex_used[j]} {sub_matrix[i,j]} # {rev_motif_dict[hex_used[i]]} x {rev_motif_dict[hex_used[j]]}\n"
                    sub_matrix_file.write(line)

    def _msa_with_characters(self, motif_dict):
        mafft_path = self.cfg.mafft_path
        chr_fasta_file = self.filename_prefix + "_motif_in_hex" 
        os.system(f"{mafft_path}/hex2maffttext {self.hex_fasta_file} > {chr_fasta_file}.ASCII")
        os.system(f"mafft --maxiterate 10 --op 1.5 --ep 0.0 --localpair --text --textmatrix {self.msa_sub_matrix}  {chr_fasta_file}.ASCII > {chr_fasta_file}_mafft_output.ASCII")
        os.system(f"{mafft_path}/maffttext2hex {chr_fasta_file}_mafft_output.ASCII > {chr_fasta_file}_mafft_output.hex")
        msa_result = list(SeqIO.parse(chr_fasta_file + "_mafft_output.hex", "fasta"))
        os.system(f"rm {chr_fasta_file}.hex")
        os.system(f"rm {chr_fasta_file}.ASCII")
        os.system(f"rm {chr_fasta_file}_mafft_output.ASCII")
        os.system(f"rm {chr_fasta_file}_mafft_output.hex")
        #os.system(f"rm {chr_fasta_file}_msa_sub_matrix.txt")



        sequence_length = int(len(msa_result[0].seq) / 2) # since each motif is represented by 2 characters
        rev_motif_chr_dict = {y: x for x, y in motif_dict.items()}
        rev_motif_chr_dict["--"] = ''

        result = numpy.zeros((len(msa_result), sequence_length),dtype=object)
        for i in range(len(msa_result)):
            single_seq = str(msa_result[i].seq)
            result[i,:] = [rev_motif_chr_dict[single_seq[i:i+2]] for i in range(0, len(single_seq), 2)]
        
        return (result, [m.description for m in msa_result])


    def _result_to_grouped_positions(self, result, msa_keys):

        nseq, sequence_length = result.shape
        unique_motifs = [set(result[:,i]) for i in range(sequence_length)]
        max_length = [0] + [max([len(e) for e in unique_motifs[i]]) for i in range(sequence_length)]
        pos = numpy.cumsum(max_length)
        start = pos[:-1]
        stop = pos[1:] 

        grouped_result = {}
        for i in range(len(msa_keys)):
            motifs = []
            for j in range(sequence_length):
                seq = result[i,j]
                if not seq:
                    continue

                if not (motifs and motifs[-1].attempt_add(seq, start=start[j])):
                    motifs.append(MotifPlacement(seq, start[j]))

            key = msa_keys[i]
            grouped_result[key] = motifs

        return grouped_result, stop[-1]


class POAAlignerNucleotide(Aligner):
    """Class for POA aligner that aligns direclty the nucleotide sequences.
       See also POAAlignerMotif for aligning the motifs directly.

    """
    def __init__(self, cfg):
        Aligner.__init__(self, cfg)


    def run(self, grouped_positions, all_unique_motifs):
        sys.stderr.write("Aligning sequences using POA Nucleotide aligner...\n"); sys.stderr.flush()
        #motif_dict, hex_used = assign_chr_to_motif(all_unique_motifs, self.printable_asii_list)
        seq_ids = []
        sequences = []
        ngrouped_positions = {}
        for sequence_id, motifs in grouped_positions.items():
            
            seq = ''.join([''.join([mp.motif] * mp.count) if not mp.singlebase else ''.join([mp.motif])   for mp in motifs])
            sequences.append(seq)
            seq_ids.append(sequence_id)
        
        aligner = pyabpoa.msa_aligner(match=1, mismatch=18, gap_open1=12, gap_ext1=4, gap_open2=36, gap_ext2=1)
        res = aligner.msa(sequences, out_cons=True, out_msa=True)
        
        result = numpy.zeros((len(seq_ids), res.msa_len),dtype='U1')
        msa_sequences = res.msa_seq[:len(seq_ids)]

        for i in range(len(seq_ids)):
            result[i,:] = list(msa_sequences[i])

        edit_distance = self._calculate_edit_distance(result, seq_ids)            

        ngrouped_positions = self._translate_position(grouped_positions, result, seq_ids)

        sequence_lengths = {key: res.msa_len for key, value in grouped_positions.items()}

        return ngrouped_positions, sequence_lengths, edit_distance


    def _translate_position(self, grouped_positions, result, seq_ids):
        """Translate the motifplacement to a new position based on a conversion dictionary from the alignment.
        
        :param grouped_positions: dict, dictionary with sequence ids as keys and list of MotifPlacement objects as values
        :param result: numpy.array, the multiple sequence alignment result. nseq * sequence_length. Each cell contains a motif or a base
        :param seq_ids: list of sequence ids

        :returns: dict, dictionary with sequence ids as keys and list of MotifPlacement objects as values
        """
        ngrouped_positions = {}
        for i in range(len(seq_ids)):
            seq_id = seq_ids[i]
            seq = result[i,:]
            keepseq = seq != '-'
            position_convert = numpy.cumsum(keepseq) - 1
            rev_pos_convert = {b:a for a,b,c in  zip(numpy.arange(len(seq)), position_convert, keepseq) if c}
            rev_pos_convert[position_convert[-1] + 1] = len(seq)

            motifs = grouped_positions[seq_id]
            nmotifs = []
            for motif in motifs:
                nmotifs.extend(motif.translate_position(rev_pos_convert, keepseq))
            ngrouped_positions[seq_id] = nmotifs
        return ngrouped_positions


class POAAlignerMotif(Aligner):
    """Class for POA aligner that aligns the motifs directly.
       Benefits:
          - shorter sequences (faster)
          - aligns always whole motifs, reducing the chance of misalignments of one motif to multiple other motifs.
       Disadvantages:
          - Does not work well in case a motif annotation needs to align to multiple single bases.
          - Cannot handle more than 23 motifs (due to mapping to the amino acids characters)
    """
    def __init__(self, cfg):
        Aligner.__init__(self, cfg)
        self.poamatrix = self.filename_prefix + "_poa_matrix" + ".txt"


    def write_poa_matrix(self, aa_letters, z, motif_dict):
        """Write the substitution matrix for the motifs to a file. 

        :param aa_letters: list of characters used in the alignment
        :param z: numpy.array, the substitution score matrix for the motifs
        :param motif_dict: dict, mapping of motifs to characters
        """

        rev_motif_dict = {y: x for x, y in motif_dict.items()}
        comment =  ' '.join([f"{aa}={rev_motif_dict[aa]}" for aa in aa_letters])
        with open(self.poamatrix, 'w') as poa_matrix_file:
            
            poa_matrix_file.write(f"# {comment}\n")
            poa_matrix_file.write(' '.join(aa_letters) + '\n')
            for i in range(len(aa_letters)):
                
                poa_matrix_file.write(' '.join([aa_letters[i]] + [str(x) for x in z[i,:]]) + '\n')

        

    def run(self, grouped_positions, all_unique_motifs):
        sys.stderr.write("Aligning sequences using POA Motif aligner...\n"); sys.stderr.flush()
        #motif_dict, hex_used = assign_chr_to_motif(all_unique_motifs, self.printable_asii_list)
        aa_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*', 'O', 'U']
        assert len(all_unique_motifs) <= len(aa_letters), "Too many motifs"
        motif_dict = dict(zip(all_unique_motifs, aa_letters))
        rev_motif_dict = {y: x for x, y in motif_dict.items()}
        aa_letters_used = list(motif_dict.values())

        z = numpy.cast[int](self._calc_sub_matrix(aa_letters_used, motif_dict))
        self.write_poa_matrix(aa_letters_used, z, motif_dict)

        
        seq_ids = []
        sequences = []
        ngrouped_positions = {}
        for sequence_id, motifs in grouped_positions.items():
            
            seq = ''.join([''.join([motif_dict[mp.motif]] * mp.count) if not mp.singlebase else ''.join([motif_dict[mc] for mc in mp.motif])   for mp in motifs])
            sequences.append(seq)
            seq_ids.append(sequence_id)
        
        aligner = pyabpoa.msa_aligner(is_aa=True, score_matrix=self.poamatrix)
        res = aligner.msa(sequences, out_cons=True, out_msa=True)
        

        result = numpy.zeros((len(seq_ids), res.msa_len),dtype=object)
        msa_sequences = res.msa_seq[:len(seq_ids)]
        rev_motif_dict['-'] = ''
        for i in range(len(seq_ids)):
            single_seq = msa_sequences[i]
            result[i,:] = [rev_motif_dict[single_seq[i]] for i in range(0, len(single_seq))]

        edit_distance = self._calculate_edit_distance(result, seq_ids)            
       
        ngrouped_positions, align_length = self._result_to_grouped_positions(result, seq_ids)

        sequence_lengths = {key: align_length for key, value in grouped_positions.items()}
        return ngrouped_positions, sequence_lengths, edit_distance
