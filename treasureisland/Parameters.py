class Parameters:
    def __init__(self):
        '''Intitialize parameters'''
        self.WINDOW_SIZE = 10000
        self.KMER_SIZE = 6
        self.UPPER_THRESHOLD = 0.80
        self.LOWER_THRESHOLD = 0.50
        self.TUNE_METRIC = 1000
        self.MINIMUM_GI_SIZE = 10000

    def set_window_size(self, window_size):
        self.WINDOW_SIZE = window_size

    def set_kmer_size(self, kmer_size):
        self.KMER_SIZE = kmer_size

    def set_upper_threshold(self, upper_threshold):
        self.UPPER_THRESHOLD = upper_threshold

    def set_lower_threshold(self, lower_threshold):
        self.LOWER_THRESHOLD = lower_threshold

    def set_tune_metric(self, tune_metric):
        self.TUNE_METRIC = tune_metric

    def set_window_size(self, minimum_gi_size):
        self.MINIMUM_GI_SIZE = minimum_gi_size

