import sys
import getopt
import gzip
from classes import ConData, file_to_con_data


def cv(argv):
    # default parameters
    impute_file_name = None
    truth_file_name = None
    
    # output parameters
    correct_string = "correct"
    wrong_string = "wrong"
    no_impute_string = "not imputed"
    no_truth_string = "imputed but no ground truth"
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "i:t:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c cv -i <impute.con> -t <truth.con> <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -i <impute.con>    imputed contacts\n")
        sys.stderr.write("  -t <truth.con>     ground truth contacts\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  leg status: correct, wrong, not imputed, imputed but no ground truth\n")
        sys.stderr.write("  tab-delimited: contact, imputed contact, ground truth contact,\n")
        sys.stderr.write("     leg 1 status, leg 2 status\n")
        return 1
    for o, a in opts:
        if o == "-i":
            impute_file_name = a
        elif o == "-t":
            truth_file_name = a
    if impute_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -i is required\n")
        return 1
    if truth_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -t is required\n")
        return 1
                       
    # read CON files
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" +  str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")

    impute_con_file = gzip.open(impute_file_name, "rb") if impute_file_name.endswith(".gz") else open(impute_file_name, "rb")
    impute_con_data = file_to_con_data(impute_con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(impute_con_data.num_cons()) + " imputed contacts (" +  str(round(100.0 * impute_con_data.num_phased_legs() / impute_con_data.num_cons() / 2, 2)) + "% legs phased)\n")

    truth_con_file = gzip.open(truth_file_name, "rb") if truth_file_name.endswith(".gz") else open(truth_file_name, "rb")
    truth_con_data = file_to_con_data(truth_con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(truth_con_data.num_cons()) + " contacts (" +  str(round(100.0 * truth_con_data.num_phased_legs() / truth_con_data.num_cons() / 2, 2)) + "% legs phased)\n")

    # sort contacts
    con_data.sort_cons()
    impute_con_data.sort_cons()
    truth_con_data.sort_cons()
    
    # analyze each contact
    for con in con_data.get_cons():
        # find corresponding contacts
        impute_cons = list(impute_con_data.get_cons_near(con, 0))
        truth_cons = list(truth_con_data.get_cons_near(con, 0))
        if len(impute_cons) != 1 or len(truth_cons) != 1:
            continue
        impute_con = impute_cons[0]
        truth_con = truth_cons[0]
        
        # determine leg 1 imputation status
        leg_1_status = None
        if con.leg_1().is_phased():
            leg_1_status = no_impute_string
        elif not truth_con.leg_1().is_phased():
            leg_1_status = no_truth_string
        elif impute_con.leg_1().get_haplotype() == truth_con.leg_1().get_haplotype():
            leg_1_status = correct_string
        else:
            leg_1_status = wrong_string

        # determine leg 2 imputation status
        leg_2_status = None
        if con.leg_2().is_phased():
            leg_2_status = no_impute_string
        elif not truth_con.leg_2().is_phased():
            leg_2_status = no_truth_string
        elif impute_con.leg_2().get_haplotype() == truth_con.leg_2().get_haplotype():
            leg_2_status = correct_string
        else:
            leg_2_status = wrong_string
                    
        sys.stdout.write("\t".join([str(con), str(impute_con), str(truth_con), leg_1_status, leg_2_status]) + "\n")
    
        
    return 0
    
