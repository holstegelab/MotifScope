import copy
import numpy as np

class MotifPlacement:
    """
        Class to store a motifplacement in a sequence. For visualization purposes.
    """
    def __init__(self, motif, start, singlebase=False):
        self.singlebase = singlebase or len(motif) == 1
        self.motif = motif
        self.start = start
        self.stop = start + len(motif)

        self.count = 1 #number of times the motif is present in the sequence. Usually int, in case of gapped alignments can be float
        self.score = None #score of the motif, used for sorting in colorbar

        self.leftborder=True #True if motifplacement starts at the start of a motif (can be False in case of gapped alignments)
        self.rightborder=True #True if motifplacement stops at the end of a motif (can be False in case of gapped alignments)

        self.first_motif_offset = 0 #if the start of motif sequence in self.motif is not at start, but before start, this is the offset (negative value) wrt. to start."

    def updateCount(self, already_done):
        """Update the count of the motif based on start/stop and motif length.
        :param already_done: float, amount of motif already present before the start position. 
                             Used if this motifplacement start halfway through a motif (i.e. in case of gapped alignment)
        """
        already_done = already_done - np.floor(already_done)
        if already_done > 0:
            self.first_motif_offset = -np.round(already_done * len(self.motif))
        else:
            self.first_motif_offset = 0

        if self.first_motif_offset != 0 and not self.singlebase:
            self.leftborder = False


        if self.singlebase:
            self.count = 1
        else:
            count = (self.stop - self.start) / len(self.motif)
            if count == np.round(count):
                self.count = int(count)
            else:
                self.count = count
        if (self.stop - (self.start + self.first_motif_offset)) % len(self.motif) != 0 and not self.singlebase:
            self.rightborder = False

    def translate_position(self, conversion_dict, keepseq):
        """Translate the motifplacement to a new position based on a conversion dictionary.
        :param conversion_dict: dict, dictionary with old positions as keys and new positions as values
        :param keepseq: np.array, boolean array with positions that in new sequence that are not gapped (due to alignment).

        Returns: list of MotifPlacement objects. Usually only one, but in case of the motifplacement crossing gaps, multiple motifplacements are created.
        

        """

        #copy the motif, adapt the start and stop positions
        nmotif = copy.deepcopy(self)
        nmotif.start = conversion_dict[self.start]
        nmotif.stop = conversion_dict[self.stop]


        if (self.stop - self.start) != (nmotif.stop - nmotif.start):
            #this is the difficult case, where the motifplacement crosses a gap or gaps
            relevant_seq = keepseq[nmotif.start:nmotif.stop]
            assert relevant_seq.sum() == self.stop - self.start
            curpos = 0
            curstart = nmotif.start
            nmotifs = []
            done_amount = 0
            while curpos < len(relevant_seq):
                if relevant_seq[curpos]:
                    curpos += 1
                else:
                    tnmotif = copy.deepcopy(nmotif)
                    tnmotif.start = curstart
                    tnmotif.stop =curpos + nmotif.start
                    #tnmotif.leftborder = False if curstart != nmotif.start and not self.singlebase else True
                    #tnmotif.rightborder = False if curpos != len(relevant_seq) and not self.singlebase else True
                    if self.singlebase:
                        tnmotif.motif = self.motif[done_amount:(done_amount + (tnmotif.stop - tnmotif.start))]
                    else:
                        tnmotif.updateCount(done_amount / len(self.motif))
                    done_amount += tnmotif.stop - tnmotif.start
                    nmotifs.append(tnmotif)

                    
                    #jump to next keep position
                    while curpos < len(relevant_seq) and not relevant_seq[curpos]:
                        curpos += 1

                    curstart = curpos + nmotif.start

            if done_amount < (self.stop - self.start):
                tnmotif = copy.deepcopy(nmotif)
                tnmotif.start = curstart
                tnmotif.stop = nmotif.stop
                #tnmotif.leftborder = False if done_amount != 0 and not self.singlebase else True
                #tnmotif.rightborder = True
                if self.singlebase:
                    tnmotif.motif = self.motif[done_amount:]
                else:
                    tnmotif.updateCount(done_amount / len(self.motif))
                nmotifs.append(tnmotif)
            #if any([not mp.leftborder or not mp.rightborder for mp in nmotifs]):
            #    import pdb; pdb.set_trace()
            return nmotifs                    
        else:
            return [nmotif]

    def attempt_add(self, motif, start=None, singlebase=False):
        """Attempt to add a motif to the current motifplacement.
        Only succeeds if the motif is the same as the current motif (or a singlebase motif), and the start position is the same as the stop position of the current motifplacement.
        
        :param motif: str, motif to add
        :param start: int, start position of the motif
        :param singlebase: bool, if True, the motifplacement is a single base motif.
        """

        if not ((self.singlebase and (singlebase or len(motif) == 1)) or self.motif == motif):
            return False

        if start is not None and self.stop != start:
            return False
       
        if self.singlebase:
            self.motif = self.motif + motif
            self.stop = self.start + len(self.motif)
        else:
            self.count += 1
            self.stop += len(motif)
        return True 

    def items(self):
        return self.start, self.stop, self.motif, self.count

    def set_score(self, score_dict):
        self.score = score_dict.get(self.motif, np.nan)

    def score_items(self):
        return self.start, self.stop, self.motif, self.count, self.score

    def __repr__(self):
        return f"({self.motif}*{self.count}:{self.start}-{self.stop})"
