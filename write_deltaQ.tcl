# Write a psf file specifying the heme atom differences in charge for Reduced and Oxidized states

package require psfgen

set heme_num 1
mol new omcs_2chains_ox_wb_ionized.psf
mol addfile omcs_2chains_ox_wb_ionized-c.pdb

mol new omcs_2chains_h${heme_num}red_wb_ionized.psf
mol addfile omcs_2chains_h${heme_num}red_wb_ionized-c.pdb
#}
#set selection_text "{resid 412 and chain Z and not name CAA HAA1 HAA2 CBA HBA1 HBA2 CGA O2A O1A CAD HAD1 HAD2 CBD HBD1 HBD2 CGD O2D O1D} or {{resid 9 12 and name SG or resid 147 13 and name CB ND1 CE1 NE2 CD2} and chain B}"
set selection_text "{resid 412 and chain Z} or {{resid 9 12 and name SG or resid 147 13 and name CB ND1 CE1 NE2 CD2} and chain B}"
#"or {resid 16 and name CB ND1 CE1 NE2 CD2 and chain B}"
set sel0 [atomselect 0 $selection_text]
set sel1 [atomselect 1 $selection_text] 
#set sel2 [atomselect 2 $selection_text]

# Set beta column values for pair interaction calculation
$sel0 set beta 1.00
$sel1 set beta 1.00
#$sel2 set beta 1.00

set ns0 [atomselect 0 "not {$selection_text}"]
set ns1 [atomselect 1 "not {$selection_text}"]
#set ns2 [atomselect 2 "not {$selection_text}"]
$ns0 set beta 2.00
$ns1 set beta 2.00
#$ns2 set beta 2.00

set indices [$sel0 get index]
set charges0 [$sel0 get charge]
set charges1 [$sel1 get charge]

set diff [vecsub $charges1 $charges0]

for {set ix 0} {$ix < [llength $indices]} {incr ix} {
	set idx [lindex $indices $ix]
	set charge [lindex $diff $ix]

	set atm1 [atomselect 0 "index $idx"]
	set atm2 [atomselect 1 "index $idx"]
	#set atm3 [atomselect 2 "index $idx"]
	$atm1 set charge $charge
	$atm2 set charge $charge
	#$atm3 set charge $charge
	#puts $charge
}

set full_system0 [atomselect 0 all]
set full_system1 [atomselect 1 all]
#set full_system2 [atomselect 2 all]

# write deltaQ psf files
$full_system0 writepsf heme_dimer/heme${heme_num}_ox_deltaQ.psf
$full_system1 writepsf heme_dimer/heme${heme_num}_red_deltaQ.psf
#$full_system2 writepsf heme_dimer/heme_pair_deltaQ.psf

# write pair interaction pdb files
$full_system0 writepdb heme_dimer/heme${heme_num}_ox_beta.pdb
$full_system1 writepdb heme_dimer/heme${heme_num}_red_beta.pdb
#$full_system2 writepdb heme_dimer/heme_pair_beta.pdb

exit
