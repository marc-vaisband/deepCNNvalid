import numpy as np
import warnings
from bamStruct import bamStruct
from alignedBAM import alignedBAM

class NoReadsError(Exception):
    """
    This is just a generic exception to signify that there are no reads at a particular position in a
    particular bam. It gets its own class so it can be caught without catching everything else.
    """
    pass


# The onehot-encoding of the nucleotide bases is coded as a global variable here.
base_enc = {
        "A": 0,
        "T": 1,
        "C": 2,
        "G": 3
    }

def onehot_encode_base(base: str):
    """
    :param base: Nucleotide as a capital letter string
    :return: np.array of length 4 which either onehot-encodes the base or is all-zero for a missing base.
    """

    result = np.zeros(4)

    try:
        index = base_enc[base]
        result[index] = 1
        return result
    except KeyError:
        warnings.warn(f"Interpreting {base} as missing base.")
        return result


class singleTensorizer:
    """
    This is in the end the class which does the heavy lifting in the tensorisation process, accessing the reads in
    the bam files, etc. It is constructed using a bamStruct object to easily pass the shape settings for the tensor.
    """


    def __init__(self, bstruct: bamStruct):
        self.window_n = bstruct.window_n
        self.window_len = 2 * self.window_n + 1
        self.max_reads = bstruct.max_reads

    def transform(self, abam: alignedBAM, position: tuple, close_file: bool = True):
        """
        This is the core transform method. It iterates through the reads and bases of the passed alignedBAM
        and onehot-encodes them before adding additional basewise information and storing everthing in a numpy array.

        :param abam: Instance of alignedBAM class containing the sequencing file to be tensorised
        :param position: A tuple of chromosome, start, stop
        :param close_file: Whether the reference to the opened aligned file should be forgotten after the transformation
        :return: Tensorised numpy.ndarray
        """

        out_data = np.zeros(shape=(self.window_len, self.max_reads, 7))

        if abam is None or abam.ID is None:
            warnings.warn("Got a None-type abam! Returning all zeros.")
            return out_data
        if not abam.is_opened:
            abam.open_file()
        if not abam.has_reads(position):
            raise NoReadsError("No reads at requested position!")


        data_mid_index = self.window_n  # This will be needed for some index juggling
        chr_name, start, stop = position

        # We iterate through all reads in the outer loop.
        for j, read in enumerate(abam.alignment_file.fetch(*position)):
            if j >= self.max_reads:
                break

            if read.is_reverse:
                read_reverse_flag = -1
            else:
                read_reverse_flag = 1

            # Index shuffling, if the read somehow does not contain the requested position, it is skipped
            read_len = len(read.seq)
            read_at_base_pos = start - read.pos - 1
            bases_to_left = start - read.pos - 1
            bases_to_right = read_len - bases_to_left
            if bases_to_left < 0 or bases_to_right < 0:
                continue

            # In the inner loop, we iterate over bases, moving outwards from the position of the substituted base.
            for m in range(self.window_n):
                r_in_data_position = data_mid_index + m
                r_in_read_position = read_at_base_pos + m
                l_in_data_position = data_mid_index - m
                l_in_read_position = read_at_base_pos - m

                # For each step, we make sure that we're still inside the read. The IndexError exceptions are to
                # not accidentally wrap around using the python -1 indexing. If the requested position index
                # exceeds the length of the read, an IndexError is also automatically raised.
                try:
                    if l_in_data_position < 0 or l_in_read_position < 0:
                        raise IndexError
                    out_data[l_in_data_position, j, :4] = onehot_encode_base(read.seq[l_in_read_position])
                    out_data[l_in_data_position, j, 4] = read.query_qualities[l_in_read_position]
                    out_data[l_in_data_position, j, 5] = read.mapping_quality
                    out_data[l_in_data_position, j, 6] = read_reverse_flag
                except IndexError:
                    pass
                try:
                    if r_in_data_position < 0 or r_in_read_position < 0:
                        raise IndexError
                    out_data[r_in_data_position, j, :4] = onehot_encode_base(read.seq[r_in_read_position])
                    out_data[r_in_data_position, j, 4] = read.query_qualities[r_in_read_position]
                    out_data[r_in_data_position, j, 5] = read.mapping_quality
                    out_data[r_in_data_position, j, 6] = read_reverse_flag
                except IndexError:
                    pass

        if close_file:
            abam.close_file()

        if not np.count_nonzero(out_data) > 0:
            raise NoReadsError("No reads at requested position!")

        return out_data
