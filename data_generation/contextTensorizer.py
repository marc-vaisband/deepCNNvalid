import numpy as np
import pandas as pd
from pairLocus import pairLocus
from singleTensorizer import singleTensorizer
from bamPhonebook import bamPhonebook


def nocontext_tensorize_plocus(stensorizer: singleTensorizer, plocus: pairLocus, abam_dict: dict):
    """
    Top-level function to obtain a no-context (i.e.: only germline and tumour are stacked) tensorisation of a candidate
    variant.

    :param stensorizer: Instance of singleTensorizer, containing all relevant tensorisation settings, which can be
    applied to alignedBAM.
    :param plocus: Instance of pairLocus, representing candidate variant
    :param abam_dict: Dictionary linking BAM IDs to alignedBAM objects
    :return: Tensorised np.ndarray with germline and tumour tracks concatenated along depth dimension
    """


    # Access pairLocus GL and CL abams and tensorize them together, stacking in depth direction
    gl_abam = abam_dict[plocus.GL_ID]
    cl_abam = abam_dict[plocus.CL_ID]
    pos = plocus.position

    return np.concatenate([stensorizer.transform(gl_abam, pos),
                           stensorizer.transform(cl_abam, pos)],
                          axis=2)


def tensorize_list(stensorizer: singleTensorizer, position: tuple, abam_list: list, stack_axis: int = 2):
    """
    Utility function to combine a list of sequencing tracks into a tensor.

    :param stensorizer: Instance of singleTensorizer, containing all relevant tensorisation settings, which can be
    applied to alignedBAM.
    :param position: Tuple of chromosome, start, stop
    :param abam_list: List of alignedBAM objects to be tensorised at given position
    :param stack_axis: Axis along which sequencing tracks are concatenated
    :return: Tensorised np.ndarray with all tracks included in abam_list concatenated along depth dimension
    """


    assert len(abam_list) > 1
    # Given a list of abams, tensorize all and stack along read direction in provided order
    return np.concatenate([stensorizer.transform(a, position) for a in abam_list],
                          axis=stack_axis)


def get_diffline_comparison_mask(metadf: pd.DataFrame, cl_line: str, cl_library: str, own_cl_ID: str):
    """
    Utility function to obtain masks of admissible comparison tracks of a different transplantation line from tumour,
    but sequenced with the same library preparation.

    :param metadf: Metadata stored in a pandas dataframe; columns including 'bam_ID', 'line_ID', 'library_ID', 'type'
    :param cl_line: ID of transplantation line of tumour
    :param cl_library: ID of library preparation of tumour
    :param own_cl_ID: ID of tumour track
    :return: Boolean mask of same length as the index of metadf, indicating tracks which are from
    """


    diffline_mask = metadf["line_ID"] != cl_line
    samelib_mask = metadf["library_ID"] == cl_library
    type_mask = metadf["type"] == "CL"
    not_self_mask = metadf["bam_ID"] != own_cl_ID

    return np.logical_and.reduce([diffline_mask, samelib_mask, type_mask, not_self_mask])


def get_sameline_comparison_mask(metadf, cl_line, cl_library, own_cl_ID):
    """
    Utility function to obtain masks of admissible comparison tracks of a same transplantation line as tumour,
    and sequenced with the same library preparation.

    :param metadf: Metadata stored in a pandas dataframe; columns including 'bam_ID', 'line_ID', 'library_ID', 'type'
    :param cl_line: ID of transplantation line of tumour
    :param cl_library: ID of library preparation of tumour
    :param own_cl_ID: ID of tumour track
    :return: Boolean mask of same length as the index of metadf, indicating tracks which are from
    """


    sameline_mask = metadf["line_ID"] == cl_line
    samelib_mask = metadf["library_ID"] == cl_library
    type_mask = metadf["type"] == "CL"
    not_self_mask = metadf["bam_ID"] != own_cl_ID

    return np.logical_and.reduce([sameline_mask, samelib_mask, type_mask, not_self_mask])


def random_context_tensorize_once(stensorizer: singleTensorizer, plocus: pairLocus, label,
                                  bambook: bamPhonebook, abam_dict: dict,
                                  k_diffline_samelib: int, k_sameline_samelib: int, stack_axis=2):
    """
    Top-level function to context-tensorize a candidate variant with a random choice of context tracks among the ones
    which are admissible.


    :param stensorizer: Instance of singleTensorizer, containing all relevant tensorisation settings, which can be
    applied to alignedBAM.
    :param plocus: Instance of pairLocus, representing candidate variant
    :param label: Class label of pairLocus, taken only for compatibility and unused
    :param bambook: Instance of bamPhonebook, storing physical locations of sequencing data for tensorisation
    :param abam_dict: abam_dict: Dictionary linking BAM IDs to alignedBAM objects
    :param k_diffline_samelib:
    :param k_sameline_samelib:
    :param stack_axis: Axis along which tracks are concatenated, should be 2 for depth concatenation as outlined in the
    manuscript.
    :return: Tensorised np.ndarray with germline, tumour and comparison tracks concatenated along stack_axis
    """

    GL_ID = plocus.GL_ID
    CL_ID = plocus.CL_ID

    cl_line = bambook.bam_metadf.loc[CL_ID, "line_ID"]
    cl_library = bambook.bam_metadf.loc[CL_ID, "library_ID"]
    try:
        assert isinstance(cl_line, str)
        assert isinstance(cl_library, str)
    except AssertionError:
        raise ValueError("Problem with parsing line/library information, they are not strings!")

    # Comparison bams from different line:
    if k_diffline_samelib > 0:
        diffline_comparison_mask = get_diffline_comparison_mask(metadf=bambook.bam_metadf, cl_line=cl_line,
                                                                cl_library=cl_library, own_cl_ID=CL_ID)
        possible_diffline_comparison_df = bambook.bam_metadf.loc[diffline_comparison_mask, :]
        for ID in possible_diffline_comparison_df.index:
            if not abam_dict[ID].has_reads(plocus.position):
                possible_diffline_comparison_df.drop(ID, axis=0, inplace=True)
    else:
        possible_diffline_comparison_df = pd.DataFrame()

    # Comparison bams from same line:
    if k_sameline_samelib > 0:
        sameline_comparison_mask = get_sameline_comparison_mask(metadf=bambook.bam_metadf, cl_line=cl_line,
                                                                cl_library=cl_library, own_cl_ID=CL_ID)
        possible_sameline_comparison_df = bambook.bam_metadf.loc[sameline_comparison_mask, :]
        for ID in possible_sameline_comparison_df.index:
            if not abam_dict[ID].has_reads(plocus.position):
                possible_sameline_comparison_df.drop(ID, axis=0, inplace=True)
    else:
        possible_sameline_comparison_df = pd.DataFrame()

    # Now make random choice for tensor, padding with 0 if not enough
    if len(possible_diffline_comparison_df) >= k_diffline_samelib:
        diffline_comp_IDs = list(
            np.random.choice(possible_diffline_comparison_df.index, size=k_diffline_samelib, replace=False))
    else:
        n_missing_diffline = k_diffline_samelib - len(possible_diffline_comparison_df)
        diffline_comp_IDs = list(possible_diffline_comparison_df.index) + [None] * n_missing_diffline

    if len(possible_sameline_comparison_df) >= k_sameline_samelib:
        sameline_comp_IDs = list(
            np.random.choice(possible_sameline_comparison_df.index, size=k_sameline_samelib, replace=False))
    else:
        n_missing_sameline = k_sameline_samelib - len(possible_sameline_comparison_df)
        sameline_comp_IDs = list(possible_sameline_comparison_df.index) + [None] * n_missing_sameline


    # Finally, tensorize and return
    out_tensor = tensorize_list(stensorizer=stensorizer, position=plocus.position, stack_axis=stack_axis,
                                abam_list=[abam_dict[ID] for ID in [GL_ID, CL_ID] +
                                           diffline_comp_IDs + sameline_comp_IDs])

    comp_list = [GL_ID, CL_ID] + diffline_comp_IDs + sameline_comp_IDs

    return out_tensor, comp_list
