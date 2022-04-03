class bamStruct:
    """
    This class holds information about the 'settings' of tensorization, namely how many tracks are stacked together,
    in which window bases around the call are considered, and how many reads per sequencing track are included.

    A composite tensor must consist at least of a GL and a CL tensor, with optional comparison tensors added
    "underneath" along the last axis.

    Instances of this class are passed over to singleTensorizer objects.
    """
    def __init__(self, window_n: int, max_reads: int, n_comparison_bams: int = 0):
        """
        Constructor

        :param window_n: Number of bases around the call start position to consider.
        :param max_reads: Maximal number of reads to consider. If fewer are present, zero-padding is applied
        :param n_comparison_bams: Number of additional sequencing tracks besides germline and tumour included.
        """


        self.window_n = window_n
        self.window_len = 2 * window_n + 1
        self.max_reads = max_reads
        self.n_comparison_bams = n_comparison_bams
        self.n_bams_total = 2 + self.n_comparison_bams
        self.single_tensor_shape = (
            self.window_len,                      # Width of "read window" considered
            max_reads * (2 + n_comparison_bams),  # 2 for GL and CL, plus all comparisons, each with max_reads reads
            7)                                    # 4 for base-onehot, then base quality, alignment quality, is_reverse
