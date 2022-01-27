class pairLocus:
    """
    This class holds the information on a locus as called by the variant caller (in our case, VarScan2).
    It is constructed from the annovar call csvs. Information about the metadata of the involved bams
    is supposed to be taken from a separate metadata structure.

    Care!!! Stop is 1-indexed for compatibility with pySAM / bamnostic.
    """
    def __init__(self, GL_ID, CL_ID, ref: str, alt: str, chromosome: str, start: int, stop: int, funcrefgene: str, gene=None, gene_detail=None, **kwargs):
        self.GL_ID = GL_ID
        self.CL_ID = CL_ID
        self.ref = ref
        self.alt = alt
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.funcrefgene = funcrefgene
        self.position = (chromosome, start, stop)
        self.gene = gene
        self.gene_detail = gene_detail

        self.idtuple = (GL_ID, CL_ID, chromosome, start, stop)
        self.kwdict = kwargs  # Probably a bad idea

        assert isinstance(start, int)
        assert isinstance(stop, int)

        if start == stop:
            raise ValueError("Stop must be at least one further than start! pySAM is 1-indexed.")

    def isSamePosition(self, pl):
        if pl.position == self.position:
            return True
        else:
            return False

    def __eq__(self, other):
        """
        Notably, this skips checking for whether the annotated category is equal!
        This allows for comparison with manually constructed ploci.
        """

        if isinstance(other, self.__class__):
            self_catless_dict = self.__dict__.copy()
            self_catless_dict.pop("category", None)

            other_catless_dict = other.__dict__.copy()
            other_catless_dict.pop("category", None)
            return self_catless_dict == other_catless_dict
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_string(self):
        return f"{self.CL_ID};{self.GL_ID};{self.chromosome};{self.start};{self.stop};{self.ref};{self.alt}"


    def write_to(self, file_pointer_object_thing, pattern=None):
        if pattern is None:
            pattern = "GL_ID;CL_ID;chromosome;start;stop;ref;alt"

        try:
            file_pointer_object_thing.write(";".join([str(self.__getattribute__(name)) for name in pattern.split(";")]) + "\n")
        except AttributeError:
            raise ValueError("Invalid attribute name specified for writing into file in pattern parameter.")




