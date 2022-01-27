import numpy as np
import warnings
from bamStruct import bamStruct

class NoReadsError(Exception):
    pass

base_enc = {
        "A": 0,
        "T": 1,
        "C": 2,
        "G": 3
    }

def onehot_encode_base(base):
    result = np.zeros(4)

    try:
        index = base_enc[base]
        result[index] = 1
        return result
    except KeyError:
        warnings.warn(f"Interpreting {base} as missing base.")
        return result


class singleTensorizer:
    def __init__(self, bstruct: bamStruct):
        self.window_n = bstruct.window_n
        self.window_len = 2 * self.window_n + 1
        self.max_reads = bstruct.max_reads

    def transform(self, abam, position: tuple, close_file = True):
        """
        abam: AlignedBAM to transform
        position: Must be a tuple of chr, start, stop, where for a SNP stop == start + 1
        """


        out_data = np.zeros(shape=(self.window_len, self.max_reads, 7))

        if abam is None or abam.ID is None:
            warnings.warn("Got a None-type abam! Returning all zeros.")
            return out_data
        if not abam.is_opened:
            abam.open_file()

        data_mid_index = self.window_n

        chr_name, start, stop = position

        assert abam.alignment_file.fetch(*position)
        i = 0
        for j, read in enumerate(abam.alignment_file.fetch(*position)):
            i += 1
            pass

        if i < 2:
            raise NoReadsError("No reads at requested position!")


        for j, read in enumerate(abam.alignment_file.fetch(*position)):
            if j >= self.max_reads:
                break

            read_len = len(read.seq)
            read_at_base_pos = start - read.pos - 1

            if read.is_reverse:
                read_reverse_flag = -1
            else:
                read_reverse_flag = 1

            bases_to_left = start - read.pos - 1
            bases_to_right = read_len - bases_to_left

            if bases_to_left < 0 or bases_to_right < 0:
                continue
            for m in np.arange(0, self.window_n):
                r_in_data_position = data_mid_index + m
                r_in_read_position = read_at_base_pos + m
                l_in_data_position = data_mid_index - m
                l_in_read_position = read_at_base_pos - m

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
                        raise ValueError("This should never ever have happened, check the single tensorizer index code")
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
