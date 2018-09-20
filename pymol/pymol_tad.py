# for TAD display
tad_chain = "X(mat)"
tad_start = 17540000
tad_end = 38080000

cmd.select("chr", "chain \"" + tad_chain + "\"")


tad_start_string = str(tad_start).rjust(9,'0')
tad_start_resn = tad_start_string[0:3]
tad_start_name = tad_start_string[3:6]
cmd.select("tad_start", "chain \"" + tad_chain + "\" and resn " + tad_start_resn + " and name " + tad_start_name)
tad_start_id = cmd.id_atom("tad_start")

tad_end_string = str(tad_end).rjust(9,'0')
tad_end_resn = tad_end_string[0:3]
tad_end_name = tad_end_string[3:6]
cmd.select("tad_end", "chain \"" + tad_chain + "\" and resn " + tad_end_resn + " and name " + tad_end_name)
tad_end_id = cmd.id_atom("tad_end")

cmd.select("tad", "id " + str(tad_start_id) + "-" + str(tad_end_id))


hide all
as lines, chr
color gray80, chr
as sticks, tad
set_bond stick_radius, 0.5, tad
spectrum b, rainbow, tad

png ~/Downloads/tad.png, 800, 600, ray=1