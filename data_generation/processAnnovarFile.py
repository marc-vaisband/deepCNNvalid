"""
This file contains utilities for parsing ANNOVAR output files (which are comma-separated).
We utilise the pandas CSV engine for convenience (to minimise dependencies, this could of course be replaced).
"""

import pandas as pd
from validatedSNPlist import build_validatedSNPlist
from pairLocus import pairLocus


def pairLocus_from_call(GL_ID: str, CL_ID: str, call: pd.Series):
    return pairLocus(GL_ID=GL_ID, CL_ID=CL_ID, ref=call["Ref"], alt=call["Alt"],
                     chromosome=call["Chr"], start=call["Start"], stop=call["End"] + 1,
                     funcrefgene=call["Func.RefGene"])

def processAnnovar(GL_ID, CL_ID, annovar_path, snplist_path, cats_used, **kwargs):
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


    return filt_ploci, labels

