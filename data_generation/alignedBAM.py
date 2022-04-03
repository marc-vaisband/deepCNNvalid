import bamnostic


class alignedBAM:
    """
    This class stores the metadata associated with a bam file and its physical path, as well as
    providing access to additional metadata.

    If speed is an issue, bamnostic (which is written entirely in Python and thus platform-independent)
    can be replaced by pySAM, however the latter may be more difficult to install, requiring certain C plugins.

    The methods open_file and close_file are implemented in case RAM is an issue. Since a single bam file might
    be invoked for many positions, it is much slower to build and destroy the AlignmentFile for every position,
    but it does mean that the dataset does not have to be held in memory.
    """

    def __init__(self, ID, bam_path, bai_path, open_immediately=False, **kwargs):
        self.bam_path = bam_path
        self.bai_path = bai_path
        self.ID = ID
        # Setting the ID and paths to None is a valid option to simulate a BAM that doesn't exist.


        # TODO: See if this is even necessary
        self.metadict = kwargs  # The kwargs can, in particular, contain the keys flowcell_ID, library_ID, ...

        self.is_opened = False
        self.alignment_file = None

        if open_immediately:
            self.open_file()

    def open_file(self):
        """
        This prompts the sequencing Alignment File to be built in memory

        :return: None
        """
        if self.ID is None:
            raise ValueError("Trying to open a None-abam!")


        if not self.is_opened:
            self.alignment_file = bamnostic.AlignmentFile(self.bam_path, index_filename=self.bai_path, mode="rb")
            self.is_opened = True
        else:
            pass


    def close_file(self):
        """
        This forgets the file in memory again, and should trigger the garbage collector to free up RAM.

        :return: None
        """
        if self.is_opened:
            self.alignment_file = None
            self.is_opened = False

    def has_reads(self, position):
        """
        This function checks whether the alignment file has any reads for a given position.
        It is very ugly, but has the advantage that it's lazy, and so doesn't need to build everything in memory.

        :param position: Tuple of chromosome, start, stop
        :return: bool
        """
        if self.ID is None:
            return False

        self.open_file()
        i = 0
        for j, read in enumerate(self.alignment_file.fetch(*position)):
            i += 1
            if i > 2:
                self.close_file()
                return True

        self.close_file()
        return False
