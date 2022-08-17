#!/bin/bash
#SBATCH --job-name=ausio2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=200MB
#SBATCH --partition=medium
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 128

time1=31250
time2=1000
time3=600000
time4=400


#################### SETUP #######################################
# lammps launch command1
LAMMPS_LAUNCH="srun /scratch/djurabek/aaleino/lammps/src/bin/lmp_g++_openmpi_trackin"
# use this 2 lines to create a dump (latter line enables appending)
DUMP_LINE="dump id all custom 20000 np_in_silica.dump id type x y z vx vy vz"
DUMP_LINE2="dump_modify id append yes"
STEPNUMBER=0
CPPSLOCATION=/scratch/djurabek/aaleino/cppscripts
#### GENERIC SCRIPTING  - DO NOT MODIFY #############


# check that mandatory files exist
for curfile in `cat temp_filedep.dat`
do
   if [ -f "$curfile" ]; then
    echo "Found $curfile."
   else 
    echo "Error: didn't find $curfile."
    exit
   fi
done
rm temp_filedep.dat

add_to_queue()
{
   noes=`ls -dv confirmation_step_*.dat -l | tail -1 | awk '{print $9}' | sed 's/[_.]/ /g' | awk '{print $3}'`
   ls > testi.dat
   echo adding. number of existing steps $noes
   rm tmp.in

   valid=${#noes}

   if [ -z "${valid}" ]
   then
      noes=0
   fi

   if [[ "${noes}" -gt "${STEPNUMBER}" ]]; then
	echo "Step  $STEPNUMBER already exists, skipping"
	STEPNUMBER=`expr $STEPNUMBER + 1`
	return
   fi


   if [ "${STEPNUMBER}" -gt "0" ]; then
   #	echo read_restart restartfile_step_${STEPNUMBER}.dat > tmp.in 
      echo read_restart restartfile.dat > tmp.in 
      echo $DUMP_LINE >> tmp.in
      echo $DUMP_LINE2 >> tmp.in
      echo 'compute		coordination all coord/atom cutoff 3.0' >> tmp.in
#      echo 'compute          ave_coord all reduce ave c_coordination' >> tmp.in
      echo 'thermo		5' >> tmp.in
      echo 'thermo_style    custom step time temp etotal press density lx ly lz dt' >> tmp.in
      echo 'thermo_modify   line one format 1 "ec  %i"' >> tmp.in
      echo 'neighbor	3.0 bin' >> tmp.in
      echo 'neigh_modify	every 10 delay 0 check yes' >> tmp.in
   fi

   cat $1 >> tmp.in
   STEPNUMBER=`expr $STEPNUMBER + 1`
   STEPNUMBERNEXT=`expr $STEPNUMBER - 1`

   echo write_restart restartfile.dat >> tmp.in 
   echo Writing backup
   echo write_restart restartfile_step_${STEPNUMBER}.dat >> tmp.in 

   echo step_number: ${STEPNUMBER} / step_name:  $1 / output in: step_${STEPNUMBER}.out / if you want to restart this step, use: restartfile_step_${STEPNUMBERNEXT}.dat >> step_names.dat

   cp tmp.in step_${STEPNUMBER}.in

   $LAMMPS_LAUNCH < tmp.in > step_${STEPNUMBER}.out

   cp step_${STEPNUMBER}.out laststep.out

# try two times again in case of fail
   if [ "$?" -ne "0" ]; then
	   $LAMMPS_LAUNCH < tmp.in > step_${STEPNUMBER}.out
	   if [ "$?" -ne "0" ]; then
		$LAMMPS_LAUNCH < tmp.in > step_${STEPNUMBER}.out 
	   	if [ "$?" -ne "0" ]; then
			echo "Program failed more than three times. Giving up. " 
			exit
		fi
	   fi
   fi

   echo done! > confirmation_step_${STEPNUMBER}.dat
}

#### END OF GENERIC SCRIPTING #############



cat > init.in << EOF

log                     ge_np_relax.log

units		        metal
atom_style	        atomic

lattice        		diamond 5.658
region        		box block 0 5 0 5 0 5

region          	mybox block -10 10.0 -10 10.0 -10 10
create_box      	3 mybox

#read silica atoms
read_dump 		       silica_au.dump 1 x y z vx vy vz box yes add yes

#define atom types
mass		                1 28.0855           # Si
mass		                2 15.999            # O
mass	                  	3 196.96657         # Au

#reset_atom_ids

# set up potentials
pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam
pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0

group		           silica type 1:2

neighbor	           0.3 bin
neigh_modify	           every 1 delay 0 check yes

# Create the inital boundary cooling zone, just so that it can be deleted later without breaking the loop.
# This will not be used!

region                     bcooling_zone block -500 500 -500 500 -4000 4000 side out units box
group                      bcool dynamic all region bcooling_zone every 100


EOF


########## SIMULATION LOGIC #####################


add_to_queue init.in

cat > keep_at_rt.in << EOF

group		           silica type 1:2
group		           gold type 3

pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam

pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0

timestep                   0.0008
velocity                   all create 1.0 6321727

fix 			   1 all nve
fix                        MYTFIX gold temp/berendsen 300.0 300.0 0.1
fix                        MYTFIX2 silica temp/berendsen 300.0 300.0 0.1
fix                        MYPFIX all press/berendsen aniso 0.0 0.0 100.0

#100 ps 

run		           $time1

dump                       createsnap all custom 100 silica_au.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes
run                        0

EOF

add_to_queue keep_at_rt.in
add_to_queue keep_at_rt.in
add_to_queue keep_at_rt.in
add_to_queue keep_at_rt.in

for impactno in `seq 1 15`
do

shift=`${CPPSLOCATION}/shiftcell_au silica_au.dump`
shiftx=`echo $shift | awk '{print $1}'`
shifty=`echo $shift | awk '{print $2}'`
    
cat > shift_atoms.in << EOF
pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam
pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0
timestep                   0.0008

fix 			   1 all nve

displace_atoms  	   all move $shiftx $shifty 0 units box

dump                       createsnap all custom 100 silica_au.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes
run                        0

EOF

add_to_queue shift_atoms.in

xmin=`$CPPSLOCATION/lammps_get_boundary_cooling_dimensions silica_au.dump 10 | awk '{print $1}'`
xmax=`$CPPSLOCATION/lammps_get_boundary_cooling_dimensions silica_au.dump 10 | awk '{print $2}'`
ymin=`$CPPSLOCATION/lammps_get_boundary_cooling_dimensions silica_au.dump 10 | awk '{print $3}'`
ymax=`$CPPSLOCATION/lammps_get_boundary_cooling_dimensions silica_au.dump 10 | awk '{print $4}'`


cat > shoot_track.in << EOF
pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam
pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0

# create the boundary cooling zone based on bash variables
# I'm not sure if these are included in the restart file, refresh just in case.

group 			   bcool delete 
region                     bcooling_zone block $xmin $xmax $ymin $ymax -4000 4000 side out units box
group                      bcool dynamic all region bcooling_zone every 100

fix                        1 all nve
fix                        bcooling bcool temp/berendsen 300.0 300.0 0.01
fix                        mytrack all trackin track.array

dump                       createsnap bcool custom 50000 border_snapshot.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes

timestep                   0.00016
run                        10

EOF

add_to_queue shoot_track.in


cat > heat_au.in << EOF
# 5 ps run to heat Au to 2000K

pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam
pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0

group 			   bcool delete 
region                     bcooling_zone block $xmin $xmax $ymin $ymax -4000 4000 side out units box
group                      bcool dynamic all region bcooling_zone every 100
group		           gold type 3

fix                        1 all nve
fix                        bcooling bcool temp/berendsen 300.0 300.0 0.01
fix 			   heat_gold gold temp/berendsen 300.0 2000.0 0.1

dump                       createsnap bcool custom 50000 border_snapshot.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes

timestep                   0.00016
run                        ${time1}

EOF

add_to_queue heat_au.in

cat > continuewithouttrack_about_100ps.in << EOF
pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam
pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0

# create the boundary cooling zone based on bash variables
# I'm not sure if these are included in the restart file, refresh just in case.

group 			   bcool delete 
region                     bcooling_zone block $xmin $xmax $ymin $ymax -4000 4000 side out units box
group                      bcool dynamic all region bcooling_zone every 100

fix                        1 all nve
fix                        bcooling bcool temp/berendsen 300.0 300.0 0.1

timestep                   0.0008
fix                        mydtfix all dt/reset 100 1.0e-5 0.001 0.02
run 			   $time3

dump                       createsnap all custom 100 silica_au.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes
run                        0

EOF

add_to_queue continuewithouttrack_about_100ps.in
cp laststep.out continuedstep.out

timesteps_left=`grep "ec " continuedstep.out | tail -1 | awk -v maxtime=$time4 '{printf("%i", (maxtime - $3) / 0.0008);}'`

initemp_au=`${CPPSLOCATION}/lammps_element_temperature silica_au.dump Au 3`
initemp_sio2=`${CPPSLOCATION}/lammps_element_temperature silica_au.dump Si 1 O 2`

cat > put_to_300K_0bar.in << EOF
pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam
pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0
timestep                   0.0008

group		           gold type 3


fix 	                   1 all nve
fix                        MYTFIX gold temp/berendsen ${initemp_au} 300.0 0.1
fix                        MYTFIX2 silica temp/berendsen ${initemp_sio2} 300.0 0.1
fix                        MYPFIX all press/berendsen aniso 0.0 0.0 100.0
run		           ${timesteps_left}

dump                       createsnap all custom 100 silica_au.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes
run                        0

EOF

add_to_queue put_to_300K_0bar.in

cat > keep_at_rt.in << EOF
pair_style                 hybrid tersoff eam morse 12.0
pair_coeff                 * * tersoff 2007_Munetoh_SiO.tersoff Si O NULL
pair_coeff                 3 3 eam Au_u3.eam

pair_coeff                 1 3 morse 0.4260    1.8139    2.5816  12.0
pair_coeff                 2 3 morse 0.0091    0.8175    5.8094  12.0

timestep                   0.0008

fix 			   1 all nve
fix                        MYTFIX gold temp/berendsen 300.0 300.0 0.1
fix                        MYTFIX2 silica temp/berendsen 300.0 300.0 0.1
fix                        MYPFIX all press/berendsen aniso 0.0 0.0 100.0

run		           $time1

dump                       createsnap all custom 100 silica_au.dump id type x y z vx vy vz
dump_modify                createsnap append no first yes
run                        0

EOF

add_to_queue keep_at_rt.in

done
