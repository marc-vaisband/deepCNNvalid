class bamStruct:
    """
    This class holds information about the order in which the tensors should be assembled.
    A composite tensor must consist at least of a GL and a CL tensor, with optional comparison tensors added
    "underneath" in the second dimension.

    It is passed to tensorizers
    """
    def __init__(self, window_n, max_reads, n_comparison_bams = 0):
        self.window_n = window_n
        self.window_len = 2 * window_n + 1
        self.max_reads = max_reads
        self.n_comparison_bams = n_comparison_bams
        self.n_bams_total = 2 + self.n_comparison_bams
        self.single_tensor_shape = (
            self.window_len,                # Width of "read window" considered
            max_reads * (2 + n_comparison_bams),  # 2 for GL and CL, plus all comparisons, each with max_reads reads
            7)                              # 4 for base onehot, then base quality, alignment quality and is_reverse
