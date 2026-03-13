import sys
import gzip


def pairs2con(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: dip-c pairs2con <contacts.pairs.gz>\n")
        sys.stderr.write("Convert hickit .pairs.gz to .con.gz format.\n")
        sys.stderr.write("Output: <input>.con.gz (replaces .pairs.gz suffix)\n")
        return 1

    input_file = argv[1]
    if input_file.endswith(".pairs.gz"):
        output_file = input_file[: -len(".pairs.gz")] + ".con.gz"
    else:
        sys.stderr.write(
            "[E::" + __name__ + "] input file must end with .pairs.gz\n"
        )
        return 1

    with gzip.open(input_file, "rt") as fin, gzip.open(
        output_file, "wt"
    ) as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            # pairs format: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1
            chr1 = fields[1]
            pos1 = fields[2]
            chr2 = fields[3]
            pos2 = fields[4]
            phase0 = fields[7] if len(fields) > 7 and fields[7] else "."
            phase1 = fields[8] if len(fields) > 8 and fields[8] else "."
            fout.write(
                "{},{},{}\t{},{},{}\n".format(chr1, pos1, phase0, chr2, pos2, phase1)
            )

    return 0
