from tqdm import tqdm


class PreprocessData:

    def __init__(self,parameters):
        self.parameters = parameters

    def generate_kmers(self, segment):
        ''' Generate overlapping kmers from a given DNA sequence
           segment : DNA sequence
        '''
        kmers = []
        for i in range(0, len(segment) - (self.parameters.KMER_SIZE - 1)):
            start = i
            end = i + self.parameters.KMER_SIZE
            kmer = segment[start:end]
            kmers.append(kmer)

        return kmers

    def get_complete_sequence_kmers(self, dna_sequence):
        '''
        Generates kmers for the entire sequence
        :param dna_sequence: DNA sequence
        :return: kmers of the entire DNA sequence
        '''
        sequence = str(dna_sequence.seq).lower()
        kmers = [self.generate_kmers(sequence)]
        return kmers

    def split_dna_sequence(self, dna_sequence):
        '''
        Divides the DNA segment into equal small segments of sizes self.window_size
        '''

        sequence = str(dna_sequence.seq).lower()
        processed_dna_seq = []
        segment_borders = []
        dna_sequence_tqdm = tqdm(range(0, len(sequence), self.parameters.WINDOW_SIZE), position=0, leave=True)
        for i in dna_sequence_tqdm:
            start = i
            if (i + self.parameters.WINDOW_SIZE) < len(sequence):
                end = i + self.parameters.WINDOW_SIZE
            else:
                end = len(sequence)
            segment = sequence[start:end]
            kmers = self.generate_kmers(segment)
            processed_dna_seq.append(kmers)
            segment_borders.append([start, end])

        return processed_dna_seq, segment_borders
