import pandas as pd
import warnings
from pairLocus import pairLocus
from validatedSNPlist import build_validatedSNPlist

"""
This file contains utilities for parsing ANNOVAR output files (which are comma-separated) into lists of pairLocus.
We utilise the pandas CSV engine for convenience (to minimise dependencies, this could of course be replaced).
"""

def pairLocus_from_call(GL_ID: str, CL_ID: str, call: pd.Series):
    """
    Constructs a pairLocus object from a single line of a pd.DataFrame representing a candidate variant.

    :param GL_ID: ID of germline sample
    :param CL_ID: ID of tumour sample
    :param call: Candidate variant represented by ANNOVAR line in a pandas Series object.
    :return: Instance of pairLocus
    """
    return pairLocus(GL_ID=GL_ID, CL_ID=CL_ID, ref=call["Ref"], alt=call["Alt"],
                     chromosome=call["Chr"], start=call["Start"], stop=call["End"] + 1,
                     funcrefgene=call["Func.RefGene"])

def processAnnovar(GL_ID: str, CL_ID: str, annovar_path, snplist_path, cats_used: list, **kwargs):
    """

    :param GL_ID: ID of germline sample
    :param CL_ID: ID of tumour sample
    :param annovar_path: Path to ANNOVAR output file of annotated candidate variants
    :param snplist_path: Path to .xlsx file of manually confirmed variants, formatted as in the documentation
                         of validatedSNPlist.build_validatedSNPlist
    :param cats_used: List of admissible entries for the 'Func.RefGene' column of ANNOVAR output.
    :param kwargs:
    :return: filt_ploci, a list of pairLocus objects representing candidate variants;
             and labels, a list of bools where 1 is a genuine mutation, and 0 is a sequencing artefact
    """

    # First, get all calls from the annovar file
    dtype_dict = {"Chr": str, "Start": int, "End": int}
    all_calls = pd.read_csv(annovar_path, sep=",", dtype=dtype_dict)

    # Construct pairLocus objects from each call, then filter for categories used
    ploci = [pairLocus_from_call(GL_ID=GL_ID, CL_ID=CL_ID, call=row) for k, row in all_calls.iterrows()]
    filt_ploci = list(filter(lambda pl: pl.funcrefgene in cats_used, ploci))

    if kwargs.get("kick_indels", False):
        filt_ploci = list(filter(lambda pl: "-" not in [pl.ref, pl.alt], filt_ploci))  # Kick indels

    # Construct list of validated mutations, get labels by associated is_present method
    vsnpl = build_validatedSNPlist(GL_ID=GL_ID, CL_ID=CL_ID, snplist_xlsx_path=snplist_path)
    labels = list(map(vsnpl.is_present, filt_ploci))
    labels = list(map(int, labels))

    return filt_ploci, labels


cats_used_global = ["exonic", "ncRNA_exonic", "splicing", "UTR3", "UTR5",
                    "UTR5;UTR3", "downstream", "upstream", "upstream;downstream"]
# This is set as a global variable because it does not change in our application, but can of course easily be made
# into a modular argument.


def get_ploci_from_annovarlist(annovarlist_path, pandas_csv_kwargs=None, **kwargs):
    """
    Top-level function which takes a csv-format list of ANNOVAR output files and their corresponding spreadsheets of
    confirmed calls, and returns tensorization-ready data.

    :param annovarlist_path: Path to csv-formatted list of ANNOVAR output files and their corresponding spreadsheets of
    confirmed calls. Must contain the columns:
            GL_bam_ID: ID of germline sample
            CL_bam_ID: ID of tumour sample
            annovar_path: path to ANNOVAR output file
            snplist_path: path to .xlsx file, formatted as in the documentation of validatedSNPlist.build_validatedSNPlist
    :param pandas_csv_kwargs: kwargs to be passed to pd.read_csv in parsing the annovarlist
    :param kwargs: Keyword arguments, passed to processAnnovar, currently only kick_indels will have an impact.
    :return: ploci, list of variant calls represented as pairLocus objects,
             labels, list of 0-1-labels where 1 is a genuine mutation, and 0 is a sequencing artefact
    """
    if pandas_csv_kwargs is None:
        pandas_csv_kwargs = {"sep": ";"}

    all_annovars_df = pd.read_csv(annovarlist_path, **pandas_csv_kwargs)
    # We utilise the pandas CSV engine for convenience

    ploci = []
    data_labels = []

    # We simply loop over all provided ANNOVAR files, and convert each into lists of pairLocus objects and labels.
    for GL_ID, CL_ID, annovar_path, snplist_path in zip(all_annovars_df["GL_bam_ID"], all_annovars_df["CL_bam_ID"],
                                                        all_annovars_df["annovar_path"], all_annovars_df["snplist_path"]):
        try:

            new_ploci, new_labels = processAnnovar(
                GL_ID=GL_ID,
                CL_ID=CL_ID,
                annovar_path=annovar_path,
                snplist_path=snplist_path,
                cats_used=cats_used_global, **kwargs)

            assert isinstance(new_labels, list)
            assert isinstance(new_ploci, list)

            data_labels += new_labels
            ploci += new_ploci
        except (ValueError, FileNotFoundError) as e:
            warnings.warn(f"Skipping the {GL_ID} vs. {CL_ID} annovar due to the following error: {e}. "
                          f"Script continues to next annovar.")
            pass

    return ploci, data_labels


def get_ploci_from_multianno(multianno_path, GL_ID: str, CL_ID: str, **kwargs):
    """
    Utility function to convert content of multianno file into variant calls represented by pairLocus instances.

    :param multianno_path: Path to .multianno file
    :param GL_ID: ID of germline sample
    :param CL_ID: ID of tumour sample
    :param kwargs: kwargs passed to constructor of pairLocus objects
    :return: List of pairLocus objects representing variant calls
    """


    with open(multianno_path, "r") as f:
        pathlines = f.read().splitlines()
    data_dict = {"Chr": [], "Start": [], "End": [], "Ref": [], "Alt": [], "Func.refGene": [], "Gene.refGene": [],
                 "GeneDetail.refGene": [], "ExonicFunc.refGene": []}

    key_order = pathlines[0].split("\t")

    for line in pathlines[1:]:
        line_parts = line.split("\t")

        for j, key in enumerate(key_order):
            if key in data_dict.keys():
                data_dict[key].append(line_parts[j])

    data_df = pd.DataFrame(data_dict)


    return [pairLocus(GL_ID=GL_ID, CL_ID=CL_ID,
                      chromosome=row["Chr"], start=int(row["Start"]), stop=int(row["End"]) + 1,
                      ref=row["Ref"], alt=row["Alt"], funcrefgene=row["Func.refGene"], gene=row["Gene.refGene"],
                      gene_detail=row["GeneDetail.refGene"], exonic_func=row["ExonicFunc.refGene"], **kwargs)
            for k, row in data_df.iterrows()]


