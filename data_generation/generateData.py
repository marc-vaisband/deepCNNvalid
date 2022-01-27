"""
This script uses provided paths to csv files tracking the available BAM files, as well as
the corresponding ANNOVAR output and lists of validated mutations. It tensorizes the
sequencing data into numpy.ndarrays and saves them as binary dumps to the stored_data folder.
"""
import numpy as np
import bamPhonebook
import alignedBAM
from bamStruct import bamStruct
from singleTensorizer import singleTensorizer
from contextTensorizer import random_context_tensorize_once, nocontext_depth_tensorize_plocus
from generatePLoci import get_ploci_from_annovarlist
import pickle
import os
import random


bambook_path = "path_here"
annovarlist_path = "path_here"

# Import the list of available bams, build alignedBAM objects for them all for later use
print("Importing list of bam files.")

bambook = bamPhonebook.build_bambook_from_csv(bambook_path)

# Now, create a dictionary for them to be accessed by ID
print("Building dict of bams.")
all_abams = {}
for bam_ID, row in bambook.bam_metadf.iterrows():
    all_abams[bam_ID] = alignedBAM.alignedBAM(ID=bam_ID, bam_path=row["bam_path"], bai_path=row["bai_path"],
                                              line_ID=row["line_ID"], library_ID=row["library_ID"])

all_abams[None] = None

print("Building ploci for SNPs...")

ploci, labels = get_ploci_from_annovarlist(annovarlist_path, kick_indels=True)
assert len(ploci) == len(labels)
print(f"Done. Length: {len(ploci)}")
with open("./data/ploci.pkl", "wb") as f:
    pickle.dump(ploci, f)


random.seed(42)
np.random.seed(42)

def tensorize_ploci(ploci_list, label_list, bamstruct, stensorizer, data_out_path, label_out_path, pl_out_path):
    pass

"""
Next, we want to build the input data. We build one dataset of only GL and CL, and
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

print("Data building done. See \"stored_data\" folder for saved files in binary format.")
