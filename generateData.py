"""
This script uses provided paths to csv files tracking the available BAM files, as well as
the corresponding ANNOVAR output and lists of validated mutations. It tensorizes the
sequencing data into numpy.ndarrays and saves them as binary dumps to the stored_data folder.
"""

import numpy as np
import pandas as pd

from data_generation.bamPhonebook import build_bambook_from_csv
from data_generation.alignedBAM import alignedBAM
from data_generation.bamStruct import bamStruct
from data_generation.singleTensorizer import singleTensorizer
from data_generation.contextTensorizer import random_context_tensorize_once, nocontext_depth_tensorize_plocus
from data_generation.generatePLoci import get_ploci_from_annovarlist, get_ploci_from_multianno
from data_generation.pairLocus import pairLocus
import pickle
import os
import random
import glob

bambook_path = "./input_data/in_facility/dummy_bamlist.csv"
annovarlist_path = "./input_data/in_facility/dummy_annovarlist.csv"

# Import the list of available bams, build alignedBAM objects for them all for later use
print("Importing list of bam files.")

bambook = build_bambook_from_csv(bambook_path)

# Now, create a dictionary for them to be accessed by ID
print("Building dict of bams.")
all_abams = {}
for bam_ID, row in bambook.bam_metadf.iterrows():
    all_abams[bam_ID] = alignedBAM(ID=bam_ID, bam_path=row["bam_path"], bai_path=row["bai_path"],
                                              line_ID=row["line_ID"], library_ID=row["library_ID"])

all_abams[None] = None

# Now that we have all bams ready, we can construct the candidate variants.
print("Building ploci for SNPs...")

ploci, labels = get_ploci_from_annovarlist(annovarlist_path, kick_indels=True)
assert len(ploci) == len(labels)
print(f"Done. Length: {len(ploci)}")
with open("./data/ploci.pkl", "wb") as f:
    pickle.dump(ploci, f)


random.seed(42)
np.random.seed(42)


"""
Next, we want to build the in-facility input data. We build one dataset of only GL and CL, and
one with up to k = 2 comparisons of same library, different line.
"""

# Tensorisation without context tracks
nocontext_folder = "../stored_data/nocontext/"
os.makedirs(os.path.abspath(nocontext_folder), exist_ok=True)

nocontext_bamstruct = bamStruct(window_n=50, max_reads=200, n_comparison_bams=0)
nocontext_stensorizer = singleTensorizer(bstruct=nocontext_bamstruct)

nocontext_depth_tensors = [nocontext_depth_tensorize_plocus(nocontext_stensorizer, pl, all_abams) for pl in ploci]
big_nocontext_tensor = np.stack(nocontext_depth_tensors, axis=0)
with open(os.path.join(nocontext_folder, "nocontext_depth_tensor.npy"), "wb") as f:
    np.save(f, big_nocontext_tensor)
print(f"Shape of big nocontext tensor is: {big_nocontext_tensor.shape}")
with open(os.path.join("nocontext_depth_labels.npy"), "wb") as f:
    np.save(f, np.array(labels).astype(int))


# For contexted data: Use a random combination of k = 2 comparisons, zeropad if not possible

contexted_folder = "../stored_data/contexted/"
os.makedirs(os.path.abspath(contexted_folder), exist_ok=True)

k2_bamstruct = bamStruct(window_n=50, max_reads=200, n_comparison_bams=2)
k2_stensorizer = singleTensorizer(bstruct=k2_bamstruct)

x1k2_tensors = []
x1k2_comp_ids = []

for pl, label in zip(ploci, labels):
    try:
        new_tensor, new_comps = random_context_tensorize_once(
            stensorizer=k2_stensorizer, stack_axis=2, plocus=pl, bambook=bambook,
            abam_dict=all_abams, k_diffline_samelib=2, k_sameline_samelib=0, label=label)
        x1k2_tensors += new_tensor
        x1k2_comp_ids += new_comps  # This is ultimately disregarded and is only tracked for debugging
    except ValueError as e:
        print(f"Tensorization in noindel part failed due to the following error: {e}. Continuing.")
        continue

big_x1k2_tensor = np.stack(x1k2_tensors, axis=0)

with open(os.path.join(contexted_folder, "x1k2_tensor.npy"), "wb") as f:
    np.save(f, big_x1k2_tensor)
print(f"Shape of big contexted tensor is: {big_x1k2_tensor.shape}")
with open(os.path.join(contexted_folder, "x1k2_labels.npy"), "wb") as f:
    np.save(f, np.array(labels).astype(int))


"""
Lastly, we also tensorise the external dataset. Due to differences in formatting, this is slightly more involved.
"""

kotani_folder = os.path.abspath("../stored_data/kotani")
os.makedirs(kotani_folder, exist_ok=True)

kotani_bamlist_path = os.path.join(kotani_folder, "bamlist.csv")
kotani_bambook = build_bambook_from_csv(kotani_bamlist_path)

kotani_abams = {}
for bam_ID, row in kotani_bambook.bam_metadf.iterrows():
    kotani_abams[bam_ID] = alignedBAM(ID=bam_ID, bam_path=row["bam_path"], bai_path=row["bai_path"],
                                              line_ID=row["line_ID"], library_ID=row["library_ID"])

kotani_abams[None] = None

kotani_ploci = []
for filepath in glob.glob("/limcr-ngs/Franz/2021-07-29_ENA_mouse_mutations_analysis/14_hc.filtered.annovar/"
                          "*.VarScan2.4.4.snp.Somatic.hc.filter.vcf.annovar.out.mm10_multianno.txt"):

    mouse_ID = filepath.split("/")[-1].split(".")[0]

    current_GL_ID = mouse_ID + "_tail"
    current_CL_ID = mouse_ID + "_BM"

    called_ploci = get_ploci_from_multianno(multianno_path=filepath,
                                            GL_ID=current_GL_ID, CL_ID=current_CL_ID, mouse_ID=mouse_ID)
    kotani_ploci += called_ploci

validated_loci_df = pd.read_csv(os.path.abspath("../input_data/kotani/kotani_snplist.csv"), sep=";")
validated_ploci_list = [pairLocus(GL_ID=row["GL_ID"], CL_ID=row["CL_ID"], ref=row["ref"], alt=row["alt"],
                                  chromosome=row["chr"], start=row["start"], stop=row["end"], funcrefgene="",
                                  gene=row["gene"]) for k, row in validated_loci_df.iterrows()]
kotani_labels = [any([vpl.isSamePosition(pl) and vpl.CL_ID == pl.CL_ID for vpl in validated_ploci_list])
                 for pl in kotani_ploci]

kotani_x1k2_tensors = []
kotani_x1k2_compIDs = []
for pl, label in zip(kotani_ploci, kotani_labels):
    try:
        new_tensor, new_comps = random_context_tensorize_once(
            stensorizer=k2_stensorizer, stack_axis=2, plocus=pl, bambook=kotani_bambook,
            abam_dict=kotani_abams, k_diffline_samelib=2, k_sameline_samelib=0, label=label)
        kotani_x1k2_tensors += new_tensor
        kotani_x1k2_compIDs += new_comps
    except ValueError as e:
        print(f"Tensorization in noindel part failed due to the following error: {e}. Continuing.")
        continue

big_kotani_x1k2_tensor = np.stack(x1k2_tensors, axis=0)

with open(os.path.join(kotani_folder, "kotani_x1k2_tensor.npy"), "wb") as f:
    np.save(f, big_kotani_x1k2_tensor)
print(f"Shape of big contexted tensor is: {big_kotani_x1k2_tensor.shape}")
with open(os.path.join(kotani_folder, "kotani_x1k2_labels.npy"), "wb") as f:
    np.save(f, np.array(kotani_labels).astype(int))

print("Data building done. See \"stored_data\" folder for saved files in binary format.")
