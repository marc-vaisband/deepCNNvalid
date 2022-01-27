import numpy as np
from processAnnovarFile import processAnnovar
from validatedSNPlist import get_manually_called_ploci
import random
import pandas as pd
import warnings
from pairLocus import pairLocus

random.seed(42)
np.random.seed(42)  # This is how reproducibility works, right? :/

"""
Strategy: First build alignedBAMs for every single BAM in the phonebook; then have a 
separate code block to mix&match them to create input tensors
"""


def get_ploci_from_annovarlist(annovarlist_path, **kwargs):

    """
    We import all the pairLoci which interest us by crawling through annovar files.
    """


    used_cats = ["exonic", "ncRNA_exonic", "splicing", "UTR3", "UTR5",
                 "UTR5;UTR3", "downstream", "upstream", "upstream;downstream"]

    all_annovars_df = pd.read_csv(annovarlist_path, sep = ";")


    ploci = []
    data_labels = []

    for GL_ID, CL_ID, annovar_path, snplist_path in zip(all_annovars_df["GL_bam_ID"], all_annovars_df["CL_bam_ID"],
                                                        all_annovars_df["annovar_path"], all_annovars_df["snplist_path"]):
        try:

            new_ploci, new_labels = processAnnovar(
                GL_ID = GL_ID,
                CL_ID = CL_ID,
                annovar_path = annovar_path,
                snplist_path = snplist_path,
                cats_used = ["exonic", "ncRNA_exonic", "splicing", "UTR3", "UTR5",
                                     "UTR5;UTR3", "downstream", "upstream", "upstream;downstream"], **kwargs)

            assert isinstance(new_labels, list)  # Just in case at some point something inexplicably becomes an array
            assert isinstance(new_ploci, list)   # and we get horrible behaviour with the +=. Yay Python!

            data_labels += new_labels
            ploci += new_ploci
        except (ValueError, FileNotFoundError) as e:
            warnings.warn(f"Skipping the {GL_ID} vs. {CL_ID} annovar due to the following error: {e}. "
                          f"Script continues to next annovar.")
            pass

    return ploci, data_labels


def get_manually_called_ploci_from_annovarlist(annovarlist_path, **kwargs):
    all_annovars_df = pd.read_csv(annovarlist_path, sep=";")

    ploci = []
    data_labels = []

    for GL_ID, CL_ID, annovar_path, snplist_path in zip(all_annovars_df["GL_bam_ID"], all_annovars_df["CL_bam_ID"],
                                                        all_annovars_df["annovar_path"],
                                                        all_annovars_df["snplist_path"]):
        try:
            new_manual_ploci = get_manually_called_ploci(GL_ID = GL_ID, CL_ID = CL_ID, snplist_xlsx_path=snplist_path)
            new_labels = [1] * len(new_manual_ploci)

            assert isinstance(new_labels, list)  # Just in case at some point something inexplicably becomes an array
            assert isinstance(new_manual_ploci, list)  # and we get horrible behaviour with the +=. Yay Python!

            data_labels += new_labels
            ploci += new_manual_ploci
        except (ValueError, FileNotFoundError) as e:
            warnings.warn(f"Skipping the {GL_ID} vs. {CL_ID} annovar due to the following error: {e}. "
                          f"Script continues to next annovar.")
            pass

    return ploci, data_labels

def get_same_position_ploci(pl: pairLocus, pl_list: list):
    return [pl_ for pl_ in pl_list if
            pl.chromosome == pl_.chromosome and pl.start == pl_.start and pl.stop == pl_.stop]