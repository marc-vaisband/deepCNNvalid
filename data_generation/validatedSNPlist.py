from pairLocus import pairLocus
import pandas as pd

class validatedSNPlist:
    """
    This class is purely for convenience of storing validated candidate variants. It is only a wrapper around
    an actual list of pairLocus objects (which store the compared samples and the position), with a slightly
    more robust __contains__ to avoid some funky behaviour.
    """
    def __init__(self, GL_ID: str, CL_ID: str, loci_list: list):
        self.GL_ID = GL_ID
        self.CL_ID = CL_ID
        self.loci_list = loci_list  # List of pairLocus objects

    def is_present(self, plocus: pairLocus):

        test_output = False
        for present_plocus in self.loci_list:
            if all([plocus.GL_ID == present_plocus.GL_ID,
                    plocus.CL_ID == present_plocus.CL_ID,
                    plocus.chromosome == present_plocus.chromosome,
                    plocus.start == present_plocus.start,
                    plocus.stop == present_plocus.stop]):

                test_output = True


        return test_output

    def __contains__(self, item: pairLocus):
        return self.is_present(item)


def build_validatedSNPlist(GL_ID, CL_ID, snplist_xlsx_path):
    """

    :param GL_ID: ID of germline sample
    :param CL_ID: ID of tumour sample
    :param snplist_xlsx_path: Path (absolute or relative) to excel file containing
    :return:
    """

    all_confirmed_calls_xlsx = pd.read_excel(snplist_xlsx_path)
    computer_confirmed_xlsx_mask = all_confirmed_calls_xlsx["comment"].apply(lambda x: x != "manually called")
    confirmed_calls_xlsx = all_confirmed_calls_xlsx[computer_confirmed_xlsx_mask]

    loci_list = []
    for k, row in confirmed_calls_xlsx.iterrows():
        loci_list.append(pairLocus(GL_ID=GL_ID, ref=row["Ref"],
                                   CL_ID=CL_ID, alt=row["Alt"],
                                   chromosome=row["Chr"],
                                   start=row["Start"],
                                   stop=row["End"] + 1,  # The +1 here is ultimately because samtools is 1-indexed
                                   funcrefgene=row["Func.RefGene"]))

    return validatedSNPlist(GL_ID=GL_ID, CL_ID=CL_ID, loci_list=loci_list)






