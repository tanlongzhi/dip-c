# chromatin near a contact
resolution = 20e3
max_sep = 10e6
hap_names = ["(pat)", "(mat)"]

leg_1 = "1,59037522,1"
leg_2 = "2,111785927,1"

leg_1 = "14,92774417,1"
leg_2 = "5,127603680,1"

leg_1 = "1,230214326,0"
leg_2 = "3,180325128,1"

leg_1 = "11,14924591,1"
leg_2 = "4,166247684,0"

leg_1 = "12,75185357,0"
leg_2 = "8,43732187,1"

leg_1 = "6,73190268,1"
leg_2 = "8,72472782,1"

chr_1, pos_1, hap_1 = leg_1.split(",")
hom_1 = chr_1 + hap_names[int(hap_1)]
pos_1 = int(pos_1)

chr_2, pos_2, hap_2 = leg_2.split(",")
hom_2 = chr_2 + hap_names[int(hap_2)]
pos_2 = int(pos_2)

start_1 = int(round((pos_1 - max_sep) / resolution) * resolution)
end_1 = int(round((pos_1 + max_sep) / resolution) * resolution)

start_2 = int(round((pos_2 - max_sep) / resolution) * resolution)
end_2 = int(round((pos_2 + max_sep) / resolution) * resolution)

cmd.select("hom_1", "chain \"" + hom_1 + "\"")
cmd.select("hom_2", "chain \"" + hom_2 + "\"")


start_1_string = str(start_1).rjust(9,'0')
cmd.select("start_1", "chain \"" + hom_1 + "\" and resn " + start_1_string[0:3] + " and name " + start_1_string[3:6])
end_1_string = str(end_1).rjust(9,'0')
cmd.select("end_1", "chain \"" + hom_1 + "\" and resn " + end_1_string[0:3] + " and name " + end_1_string[3:6])
cmd.select("subchain_1", "id " + str(cmd.id_atom("start_1")) + "-" + str(cmd.id_atom("end_1")))

start_2_string = str(start_2).rjust(9,'0')
cmd.select("start_2", "chain \"" + hom_2 + "\" and resn " + start_2_string[0:3] + " and name " + start_2_string[3:6])
end_2_string = str(end_2).rjust(9,'0')
cmd.select("end_2", "chain \"" + hom_2 + "\" and resn " + end_2_string[0:3] + " and name " + end_2_string[3:6])
cmd.select("subchain_2", "id " + str(cmd.id_atom("start_2")) + "-" + str(cmd.id_atom("end_2")))

hide all
color blue, hom_1
color red, hom_2
#as lines, hom_1 + hom_2
as sticks, subchain_1 + subchain_2
#spectrumany count, white blue, subchain_1
#spectrumany count, red black, subchain_2
set_bond stick_radius, 0.5, subchain_1 + subchain_2


hide all
#as lines, hom_1
as sticks, subchain_1
png ~/Downloads/chain_1.png, 600, 600, ray=1
#as lines, hom_2
as sticks, subchain_2
png ~/Downloads/chains.png, 600, 600, ray=1


print hom_1 + ": " + str(int(round(start_1/1e6))) + " – " + str(int(round(end_1/1e6))) + " Mb"
print hom_2 + ": " + str(int(round(start_2/1e6))) + " – " + str(int(round(end_2/1e6))) + " Mb"
