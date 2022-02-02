import glob
import numpy as np
import pandas as pd
import os
import sys
sys.path.append(os.path.abspath("../data_generation/"))
from bamStruct import bamStruct
from alignedBAM import alignedBAM
from singleTensorizer import singleTensorizer
from contextTensorizer import random_context_tensorize
from processAnnovarFile import pairLocus_from_call
from bamPhonebook import build_bambook_from_csv
from process_tapsep_multianno import pl_list_from_multianno
import pickle


"""
Build bamlist
"""

kotani_folder = os.path.abspath("./stored_data/kotani")
os.makedirs(kotani_folder, exist_ok=True)

kotani_bamlist_path = os.path.join(kotani_folder, "bamlist.csv")
with open(kotani_bamlist_path, "w") as f:
    f.write("bam_ID;line_ID;flowcell_ID;library_ID;bam_path;bai_path\n")


    for filepath in glob.glob("input_data/kotani/*.calmd.bam"):
        filename = filepath.split("/")[-1]
        associated_bai = filepath.replace(".bam", ".bai")

        bam_ID = filename.split(".")[0]
        line_ID = str(bam_ID)[0]
        flowcell_ID = "NA"
        library_ID = 1

        f.write(f"{bam_ID};{line_ID};{flowcell_ID};{library_ID};{filepath};{associated_bai}\n")

comparison_bambook = build_bambook_from_csv(kotani_bamlist_path)


"""
Build list of called ploci
"""
all_called_ploci = []
labels = []

for filepath in glob.glob("/limcr-ngs/Franz/2021-07-29_ENA_mouse_mutations_analysis/14_hc.filtered.annovar/"
                          "*.VarScan2.4.4.snp.Somatic.hc.filter.vcf.annovar.out.mm10_multianno.txt"):

    mouse_ID = filepath.split("/")[-1].split(".")[0]

    current_GL_ID = mouse_ID + "_tail"
    current_CL_ID = mouse_ID + "_BM"

    called_ploci = pl_list_from_multianno(multianno_path=filepath, GL_ID=current_GL_ID, CL_ID=current_CL_ID, mouse_ID=mouse_ID)
    all_called_ploci += called_ploci

print(f"Total number of called ploci:  {len(all_called_ploci)}")

exonic_ploci = [pl for pl in all_called_ploci if pl.funcrefgene == "exonic"]



print(f"Of them exonic: {len(exonic_ploci)}")

bams_mentioned = []
for pl in all_called_ploci:
    bams_mentioned += [pl.GL_ID, pl.CL_ID]

bams_mentioned = set(bams_mentioned)


for bam_ID in bams_mentioned:
    try:
        _ = comparison_bambook.bam_metadf.loc[bam_ID, :]

    except KeyError:
        print(f"No bam found in bamlist for ID {bam_ID}")

with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/exonic_ploci_noindels_testset1.pkl", "wb") as f:
    pickle.dump(exonic_ploci, f)

# with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/exonic_ploci_noindels_testset1.pkl", "rb") as f:
#     exonic_ploci = pickle.load(f)


"""
Build labels for ploci
"""
dtype_dict = {"Start": int, "End": int, "Chromosome": str}
all_confirmed_snps = pd.read_csv("/limcr-ngs/Franz/2021-07-29_ENA_mouse_mutations_analysis/Supp table 3_SNP list_2019_Kotani.txt", skiprows=1, sep=" ", dtype=dtype_dict)

labels = []

for pl in exonic_ploci:

    inspected_df = all_confirmed_snps[all_confirmed_snps["Mouse"] + "_tail" == pl.GL_ID]
    inspected_df = inspected_df[inspected_df["Chromosome"] == pl.chromosome]

    if pl.start in list(inspected_df["Start"]):
        labels.append(1)
    else:
        labels.append(0)


with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/exonic_labels_noindels_testset1.pkl", "wb") as f:
    pickle.dump(labels, f)

print(f"Total mutations: {sum(labels)}, artifacts: {len(labels) - sum(labels)}")


kick_regions = [("chr9", 44803354, 44881274),
                ("chr4", 87769924, 87791965)
                ]

keep = [not any([pl.chromosome == chrom and pl.start >= start and pl.stop <= stop for chrom, start, stop in kick_regions])

        for pl in exonic_ploci]

print(f"Number of ploci kicked: {len(exonic_ploci) - sum(keep)}")
subset_ploci = [pl for pl, flag in zip(exonic_ploci, keep) if flag]
with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/subset_ploci.pkl", "wb") as f:
    pickle.dump(subset_ploci, f)


keep_nonsyn = [pl.kwdict["exonic_func"] != "synonymous SNV" for pl in subset_ploci]


nonsyn_ploci = [pl for pl, flag in zip(subset_ploci, keep_nonsyn) if flag]

comparison_data_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/depth_tensor.npy"
comparison_label_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/labels.npy"
with open(comparison_data_path, "rb") as f:
    valid_data = np.load(f)
with open(comparison_label_path, "rb") as f:
    valid_labels = np.load(f)
assert len(exonic_ploci) == len(valid_data)
subset_valid_data = valid_data[keep]
subset_valid_labels = valid_labels[keep]


assert len(subset_valid_data) == len(subset_ploci)

subset_comparison_data_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/subset_depth_tensor.npy"
subset_comparison_label_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/subset_labels.npy"
with open(subset_comparison_data_path, "wb") as f:
    np.save(f, subset_valid_data)
with open(subset_comparison_label_path, "wb") as f:
    np.save(f, subset_valid_labels)


nonsyn_data = subset_valid_data[keep_nonsyn]
nonsyn_labels = subset_valid_labels[keep_nonsyn]

print(f"Number of ploci kicked for being synonymous: {len(subset_ploci) - len(nonsyn_ploci)}")

nonsyn_comparison_data_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/nonsyn_depth_tensor.npy"
nonsyn_comparison_label_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/nonsyn_labels.npy"
with open(nonsyn_comparison_data_path, "wb") as f:
    np.save(f, nonsyn_data)
with open(nonsyn_comparison_label_path, "wb") as f:
    np.save(f, nonsyn_labels)

with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/nonsyn_ploci.pkl", "wb") as f:
    pickle.dump(nonsyn_ploci, f)



#
# print("Building abams.")
# all_abams = {}
# for bam_ID, row in comparison_bambook.bam_metadf.iterrows():
#     all_abams[bam_ID] = alignedBAM(ID = bam_ID, bam_path = row["bam_path"], bai_path = row["bai_path"],
#                                               line_ID = row["line_ID"], library_ID = row["library_ID"])
#
#
# k2_bamstruct = bamStruct(window_n = 50, max_reads=200, n_comparison_bams=2)
# k2_stensorizer = singleTensorizer(bstruct=k2_bamstruct)
# x1k2_noindels_depth_tensors = []
#
# for pl, label in zip(exonic_ploci, labels):
#     try:
#         new_tensor, _, new_comps, __ = random_context_tensorize(stensorizer=k2_stensorizer, axis=2, plocus=pl, bambook=comparison_bambook, abam_dict=all_abams, k_bams=2, n_points=1, label=label)
#         assert isinstance(new_tensor, list)
#         x1k2_noindels_depth_tensors += new_tensor
#     except (ValueError) as e:
#         print(f"Tensorization in noindel part failed due to the following error: {e}. Continuing.")
#         raise
#
# x1k2_noindels_depth_tensor = np.stack(x1k2_noindels_depth_tensors, axis=0)
#
# testset_file_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1"
# os.makedirs(testset_file_path, exist_ok=True)
#
# with open(os.path.join(testset_file_path, "depth_tensor.npy"), "wb") as f:
#     np.save(f, x1k2_noindels_depth_tensor)
# with open(os.path.join(testset_file_path, "labels.npy"), "wb") as f:
#     np.save(f, np.array(labels).astype(int))


# This is a comment






