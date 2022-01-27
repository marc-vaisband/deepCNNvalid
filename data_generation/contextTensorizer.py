import numpy as np
from itertools import combinations
from pairLocus import pairLocus
from singleTensorizer import singleTensorizer, NoReadsError
from bamPhonebook import bamPhonebook


def nocontext_tensorize_plocus(stensorizer: singleTensorizer, plocus: pairLocus, abam_dict: dict):
    # Access pairLocus GL and CL abams and tensorize them together, stacking in read direction
    gl_abam = abam_dict[plocus.GL_ID]
    cl_abam = abam_dict[plocus.CL_ID]
    pos = plocus.position

    return np.concatenate([stensorizer.transform(gl_abam, pos),
                           stensorizer.transform(cl_abam, pos)],
                          axis=1)


def nocontext_depth_tensorize_plocus(stensorizer: singleTensorizer, plocus: pairLocus, abam_dict: dict):
    # Access pairLocus GL and CL abams and tensorize them together, stacking in read direction
    gl_abam = abam_dict[plocus.GL_ID]
    cl_abam = abam_dict[plocus.CL_ID]
    pos = plocus.position

    return np.concatenate([stensorizer.transform(gl_abam, pos),
                           stensorizer.transform(cl_abam, pos)],
                          axis=2)


def tensorize_list(stensorizer: singleTensorizer, position: tuple, abam_list: list, axis: int = 1):
    assert len(abam_list) > 1  # Otherwise, should check what happens when concatenate is called on just one array
    # Given a list of abams, tensorize all and stack along read direction in provided order
    return np.concatenate([stensorizer.transform(a, position) for a in abam_list],
                          axis=axis)


def simple_context_tensorize(stensorizer: singleTensorizer, plocus: pairLocus,
                             bambook: bamPhonebook, abam_dict: dict, k: int, axis=1):
    # Try to find k comparisons with same library, different line, if not possible, pad with None, then tensorize
    # Returns one data tensor!

    GL_ID = plocus.GL_ID
    CL_ID = plocus.CL_ID
    try:
        cl_line = bambook.bam_metadf.loc[CL_ID, "line_ID"]
        cl_library = bambook.bam_metadf.loc[CL_ID, "library_ID"]
    except KeyError:
        raise KeyError(f"Unknown line or library for sample {CL_ID}")

    # First step: See how many possible comparisons there are
    possible_comparison_mask = np.logical_and(np.array(bambook.bam_metadf["line_ID"]) != cl_line,
                                              np.array(bambook.bam_metadf["library_ID"]) == cl_library)

    possible_comparison_df = bambook.bam_metadf.loc[possible_comparison_mask, :]

    for ID in possible_comparison_df.index:
        try:
            _ = stensorizer.transform(abam_dict[ID], plocus.position)
        except NoReadsError:
            possible_comparison_df.drop(ID, axis=0, inplace=True)

    if (n_possible_comparisons := len(possible_comparison_df)) < k:
        # Take whatever is there, then pad to k with zeros/Nones
        n_pads = k - n_possible_comparisons

        comparison_IDs = list(possible_comparison_df.index) + [None] * n_pads
        print(comparison_IDs)

        return tensorize_list(stensorizer=stensorizer, position=plocus.position,
                              abam_list=[abam_dict[ID] for ID in [GL_ID, CL_ID] +
                                         possible_comparison_df.index] + [None] * n_pads), comparison_IDs

    else:
        comparison_IDs = np.random.choice(possible_comparison_df.index, size=k, replace=False)
        return tensorize_list(stensorizer=stensorizer, position=plocus.position,
                              abam_list=[abam_dict[ID] for ID in [GL_ID, CL_ID] + comparison_IDs],
                              axis=axis), comparison_IDs


def get_diffline_comparison_mask(metadf, cl_line, cl_library, own_cl_ID):
    diffline_mask = metadf["line_ID"] != cl_line
    samelib_mask = metadf["library_ID"] == cl_library
    type_mask = metadf["type"] == "CL"
    not_self_mask = metadf["bam_ID"] != own_cl_ID

    return np.logical_and.reduce([diffline_mask, samelib_mask, type_mask, not_self_mask])


def get_sameline_comparison_mask(metadf, cl_line, cl_library, own_cl_ID):
    sameline_mask = metadf["line_ID"] == cl_line
    samelib_mask = metadf["library_ID"] == cl_library
    type_mask = metadf["type"] == "CL"
    not_self_mask = metadf["bam_ID"] != own_cl_ID

    return np.logical_and.reduce([sameline_mask, samelib_mask, type_mask, not_self_mask])


def random_context_tensorize_once(stensorizer: singleTensorizer, plocus: pairLocus, label,
                                  bambook: bamPhonebook, abam_dict: dict, k_diffline_samelib: int, k_sameline_samelib: int,
                                  stack_axis=2):
    GL_ID = plocus.GL_ID
    CL_ID = plocus.CL_ID

    cl_line = bambook.bam_metadf.loc[CL_ID, "line_ID"]
    cl_library = bambook.bam_metadf.loc[CL_ID, "library_ID"]
    try:
        assert isinstance(cl_line, str)
        assert isinstance(cl_library, str)
    except AssertionError:
        print(cl_line, cl_library)
        print(type(cl_line), type(cl_library))
        raise

    # Comparison bams from different line:
    diffline_comparison_mask = get_diffline_comparison_mask(metadf=bambook.bam_metadf, cl_line=cl_line,
                                                            cl_library=cl_library, own_cl_ID=CL_ID)
    possible_diffline_comparison_df = bambook.bam_metadf.loc[diffline_comparison_mask, :]
    for ID in possible_diffline_comparison_df.index:
        if not abam_dict[ID].has_reads(plocus.position):
            possible_diffline_comparison_df.drop(ID, axis=0, inplace=True)

    # Comparison bams from same line:
    sameline_comparison_mask = get_sameline_comparison_mask(metadf=bambook.bam_metadf, cl_line=cl_line,
                                                            cl_library=cl_library, own_cl_ID=CL_ID)
    possible_sameline_comparison_df = bambook.bam_metadf.loc[sameline_comparison_mask, :]
    for ID in possible_sameline_comparison_df.index:
        if not abam_dict[ID].has_reads(plocus.position):
            possible_sameline_comparison_df.drop(ID, axis=0, inplace=True)

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

    out_tensor = tensorize_list(stensorizer=stensorizer, position=plocus.position, axis=stack_axis,
                                abam_list=[abam_dict[ID] for ID in [GL_ID, CL_ID] +
                                           diffline_comp_IDs + sameline_comp_IDs])

    comp_list = [GL_ID, CL_ID] + diffline_comp_IDs + sameline_comp_IDs

    return out_tensor, comp_list
