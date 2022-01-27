import pandas as pd
from bamStruct import bamStruct


class bamPhonebook:
    """
    This is basically just a pandas dataframe with some utility functions stapled on top in order to
    find IDs from paths, paths from IDs, and so on.
    """
    def __init__(self, bam_metadf: pd.DataFrame):
        self.bam_metadf = bam_metadf
        self.bam_metadf["bam_ID"] = self.bam_metadf.index
        self.bam_metadf["type"] = self.bam_metadf["bam_ID"].apply(lambda ID: "CL" if ID[0] == "C" else "GL")

    def find_path(self, bam_ID):
        return self.bam_metadf.loc[bam_ID, "path"]


def build_bambook_from_csv(path, sep=";", index_col=0, **kwargs):
    return bamPhonebook(pd.read_csv(path, sep=sep, index_col=index_col, **kwargs,
                                    dtype={"bam_ID": str, "library_ID": str, "line_ID": str}))
