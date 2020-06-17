#!/usr/bin/env bash

for filename in ${PWD}/parameter_values/*.csv; do
    # get number of lines (parameter value combinations) in file
    num_lines=$(wc -l < "$filename")
    for ((param_combo=1; param_combo<=num_lines; param_combo++)); do
	line=$(sed -n ''"$param_combo"'p' $filename | tr -d "\"")
	IFS=',' read -ra PARAM_ARRAY <<< "$line"
	
	case $filename in
	    *HSE_*.csv)
		# 4 parameters (population_size, number_generations, number_replicates, selection_coefficient)
		if [ "${#PARAM_ARRAY[@]}" == "4" ]; then
		    ps=${PARAM_ARRAY[0]}
		    ng=${PARAM_ARRAY[1]}
		    nr=${PARAM_ARRAY[2]}
		    sc=${PARAM_ARRAY[3]}
		    # run python script
		    cd scripts/
		    python3 run_haploid_single_environment.py $ps $ng $nr $sc
		    cd ..
		else
		    echo "incorrect number of parameters in $filename"
		fi
		;;
	    *DSE_*.csv)
		# 5 parameters (population_size, number_generations, number_replicates,
		# selection_coefficient_homozygote, selection_coefficient_heterozygote)
		if [ "${#PARAM_ARRAY[@]}" == "5" ]; then
		    ps=${PARAM_ARRAY[0]}
		    ng=${PARAM_ARRAY[1]}
		    nr=${PARAM_ARRAY[2]}
		    scho=${PARAM_ARRAY[3]}
		    sche=${PARAM_ARRAY[4]}
		    # run python script
		    cd scripts/
		    python3 run_diploid_single_environment.py $ps $ng $nr $scho $sche
		    cd ..
		else
		    echo "incorrect number of parameters in $filename"
		fi
		;;
	    *HTE_*.csv)
		# 8 parameters (population_size, gen_env_1, gen_env_2, number_replicates,
		# selection_coefficient_A_env_1, selection_coefficient_A_env_2,
		# selection_coefficient_a_env_1, selection_coefficient_a_env_2)
		if [ "${#PARAM_ARRAY[@]}" == "8" ]; then
		    ps=${PARAM_ARRAY[0]}
		    ge1=${PARAM_ARRAY[1]}
		    ge2=${PARAM_ARRAY[2]}
		    nr=${PARAM_ARRAY[3]}
		    scAe1=${PARAM_ARRAY[4]}
		    scAe2=${PARAM_ARRAY[5]}
		    scae1=${PARAM_ARRAY[6]}
		    scae2=${PARAM_ARRAY[7]}
		    # run python script
		    cd scripts/
		    python3 run_haploid_two_environments.py $ps $ge1 $ge2 $nr $scAe1 $scAe2 $scae1 $scae2
		    cd ..
		else
		    echo "incorrect number of parameters in $filename"
		fi
		;;
	    *HTEOE_*.csv)
		# 7 parameters (population_size, gen_env_1, gen_env_2, number_replicates,
		# selection_coefficient_A_env_1, selection_coefficient_A_env_2,
		# selection_coefficient_a_env_1, selection_coefficient_a_env_2)
		if [ "${#PARAM_ARRAY[@]}" == "7" ]; then
		    ps=${PARAM_ARRAY[0]}
		    ng=${PARAM_ARRAY[1]}
		    nr=${PARAM_ARRAY[2]}
		    scA1=${PARAM_ARRAY[3]}
		    scA2=${PARAM_ARRAY[4]}
		    sca1=${PARAM_ARRAY[5]}
		    sca2=${PARAM_ARRAY[6]}
		    # run python script
		    cd scripts/
		    python3 run_haploid_two_effects_one_environment.py $ps $ng $nr $scA1 $scA2 $sca1 $sca2
		    cd ..
		else
		    echo "incorrect number of parameters in $filename"
		fi
		;;
	    *)
		echo "should not execute, error with case switching"
		;;
	esac
	
    done
done

# analyse output
cd ../../analysis/metric_calculation/
matlab -nodisplay -r "run process_persistence_probs_to_metric.m; quit"
cd ../../run/population_genetics_models/
