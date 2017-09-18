import sys
import getopt
import gzip
import copy
from classes import ConData, file_to_con_data, Leg, Par, ParData

def impute(argv):
    # default parameters
    num_rounds = 3
    max_impute_distance = 10000000
    min_impute_votes = 3
    min_impute_vote_fraction = 0.9
    max_intra_hom_separation = 10000000
    min_inter_hom_separation = 100000000
    max_clean_distance = 10000000
    min_clean_count = 2
    is_male = False
    par_data = None
    
    # presets
    h_par = ParData("X", "Y")
    h_par.add_par(Par("X", 60000, 2699520, "Y", 10000))
    h_par.add_par(Par("X", 154931043, 155260560, "Y", 59034049))

    m_par = ParData("chrX", "chrY")
    m_par.add_par(Par("chrX", 169969758, 170931299, "chrY", 90745844))
    
    presets = {
        "f":[False, None],
        "hm":[True, h_par],
        "mm":[True, m_par]}
    preset_descriptions = {
        "f":"female",
        "hm":"human male (hg19, no \"chr\" prefix)",
        "mm":"mouse male (mm10)"}
    
    # progress display parameters
    display_max_num_legs = 20
    display_num_cons = 10000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "d:v:f:D:C:p:r:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac impute [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -d INT     max distance (bp, L-1/2 norm) for imputing haplotypes [" + str(max_impute_distance) + "]\n")
        sys.stderr.write("  -v INT     min votes for the major haplotypes [" + str(min_impute_votes) + "]\n")
        sys.stderr.write("  -f FLOAT   min vote fraction for the major haplotypes [" + str(min_impute_vote_fraction) + "]\n\n")
        sys.stderr.write("  -s INT     max separation (bp) for assuming intra-homologous contacts [" + str(max_intra_hom_separation) + "]\n")
        sys.stderr.write("  -S INT     min separation (bp) for allowing inter-homologous contacts [" + str(min_inter_hom_separation) + "]\n\n")
        sys.stderr.write("  -r INT     num of rounds of imputing partially phased contacts [" + str(num_rounds) + "]\n\n")
        sys.stderr.write("  -D INT     max distance (bp, L-1/2 norm) for removing isolated contacts [" + str(max_clean_distance) + "]\n")
        sys.stderr.write("  -C INT     min neighbor count for an unisolated contact [" + str(min_clean_count) + "]\n\n")
        sys.stderr.write("  -p STR     presets for PARs and sex: [f]\n")
        for preset in sorted(presets.keys()):
            sys.stderr.write("               " + preset + " = " + preset_descriptions[preset] + "\n")
        return 1
    for o, a in opts:
        if o == "-d":
            max_impute_distance = int(a)
        elif o == "-v":
            min_impute_votes = int(a)
        elif o == "-f":
            min_impute_vote_fraction = float(a)
        elif o == "-r":
            num_rounds = int(a)
        elif o == "-D":
            max_clean_distance = int(a)
        elif o == "-C":
            min_clean_count = int(a)  
        elif o == "-p":
            try:
                is_male, par_data = presets[a]
                sys.stderr.write("[M::" + __name__ + "] use preset " + a + " = " + preset_descriptions[a] + "\n")
            except KeyError:
                sys.stderr.write("[E::" + __name__ + "] unknown preset\n")
                return 1   
                            
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    
    if is_male:
        # discard all contacts in PARs
        con_data.clean_in_par(par_data)
        # set haplotypes for all the remaining X (mat) or Y (pat) contacts
        con_data.set_non_par_hap_tuple_male(par_data)
    sys.stderr.write("[M::" + __name__ + "] removed PARs for male, leaving " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
            
    # find two subsets: A. phased or partly phased; B. completely unphased, inter-chromosomal
    no_phased_inter_chr_con_data = ConData()
    some_phased_con_data = ConData()
    for con in con_data.get_cons():
        if con.num_phased_legs() > 0:
            some_phased_con_data.add_con(con)
        elif not con.is_intra_chr():
            no_phased_inter_chr_con_data.add_con(con)
    sys.stderr.write("[M::" + __name__ + "] found " + str(no_phased_inter_chr_con_data.num_cons()) + " completely unphased, inter-chromosomal contacts\n")
    sys.stderr.write("[M::" + __name__ + "] found " + str(some_phased_con_data.num_cons()) + " phased or partly phased contacts\n")
            
    # pass 1: impute A with A
    some_phased_con_data.sort_cons()
    evidence_con_data = copy.deepcopy(some_phased_con_data)
    for round_id in range(num_rounds):
        target_con_data = copy.deepcopy(some_phased_con_data) # always use the original as the target of imputation
        sys.stderr.write("[M::" + __name__ + "] pass 1 round " + str(round_id + 1) + ": imputing partly phased contacts\n")
        
        # impute
        evidence_con_data.sort_cons()
        target_con_data.impute_from_con_data(evidence_con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation, min_inter_hom_separation)
        
        sys.stderr.write("[M::" + __name__ + "] pass 1 round " + str(round_id + 1) + " done: imputed " + str(target_con_data.num_phased_cons()) + " contacts (" + str(round(100.0 * target_con_data.num_phased_cons() / target_con_data.num_cons(), 2)) + "% of all phased or partly phased contacts)\n")
        evidence_con_data = target_con_data # use the result of imputation to refine the original in the next round
        
    some_phased_con_data = target_con_data
    
    # pass 2: clean imputed A
    some_phased_con_data.sort_cons()
    some_phased_con_data.clean_unphased()
    pass_1_num_cons = some_phased_con_data.num_cons()
    some_phased_con_data.clean_isolated_phased(copy.deepcopy(some_phased_con_data), max_clean_distance, min_clean_count)
    pass_2_num_cons = some_phased_con_data.num_cons()
    sys.stderr.write("[M::" + __name__ + "] pass 2 done: removed " + str(pass_1_num_cons - pass_2_num_cons) + " isolated contacts (" + str(round(100.0 * (pass_1_num_cons - pass_2_num_cons) / pass_1_num_cons, 2)) + "%)\n")

    # pass 3: impute B with A
    no_phased_inter_chr_con_data.impute_from_con_data(some_phased_con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation, min_inter_hom_separation)
    sys.stderr.write("[M::" + __name__ + "] pass 3 done: imputed " + str(no_phased_inter_chr_con_data.num_phased_cons()) + " contacts (" + str(round(100.0 * no_phased_inter_chr_con_data.num_phased_cons() / no_phased_inter_chr_con_data.num_cons(), 2)) + "% of all completely unphased, inter-chromosomal contacts)\n")
    no_phased_inter_chr_con_data.clean_unphased()
    some_phased_con_data.merge_with(no_phased_inter_chr_con_data)
    
    # write output
    con_data = some_phased_con_data
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    sys.stdout.write(con_data.to_string()+"\n")

    
    return 0
    