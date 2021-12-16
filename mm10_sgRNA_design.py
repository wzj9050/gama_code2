import argparse
parser = argparse.ArgumentParser()
parser.add_argument('k',type=int,help='the length of sgRNA')
parser.add_argument('--gap',type=int,default=3,
                    help='the distance between the cleavage position and the start of PAM')
parser.add_argument('--compensation',type=int,default=1,
                    help='compensate the N in NGG(for more details please read the introduction)')
parser.add_argument('gtf_path',type=str,help='the path of reference annotations')
parser.add_argument('fasta_path',type=str,help='the path of reference gene file in fasta')
parser.add_argument('-cds','--cdf_scoring_file',type=str,default=r'..\res\STable 19 FractionActive_dlfc_lookup(1)(1).xlsx',
                    help='the path of cds scoring criteria sotred in .xlsx')
parser.add_argument('-f','--filter_threshold',type=int,default=4,
                    help='the threshold to filter these sequences share the threshold length bases ahead to increase alignment')
parser.add_argument('-ali_dict','--alignment_scoring_dictionary',type=dict,
                    default={"match": 3, "mismatch": -2, 'gap_opening': -10,'gap_extension': -1},
                    help='the scoring criteria used in the local alignment scoring matrix')
