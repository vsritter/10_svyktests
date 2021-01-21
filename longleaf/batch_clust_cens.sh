
# Variables used on the R script
module load r/4.0.3

array_env_simid=(1 2)
array_env_seed=(68180 57782)
array_env_s_size=(500 500)
array_env_cens=(20 20)
array_env_psu_size=(5 5)
array_env_effect=(0 1)

env_test=(pepe gray)
env_rep=5

export ENV_REP=${env_rep}
for t in ${!env_test[@]}; do
  export ENV_TEST=${env_test[t]}
  for idx in ${!array_env_simid[@]}; do
    export ENV_SIMID=${array_env_simid[idx]}
    export ENV_SEED=${array_env_seed[idx]}
    export ENV_S_SIZE=${array_env_s_size[idx]}
    export ENV_CENS=${array_env_cens[idx]}
    export ENV_PSU_SIZE=${array_env_psu_size[idx]}
    export ENV_EFFECT=${array_env_effect[idx]}
    sbatch call_sim.sh
  done
done














