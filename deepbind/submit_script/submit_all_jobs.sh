set -o nounset -o pipefail -o errexit

IFS=$'\n'
for model in $(cat ../data/test_model.ids)
do
	# reformat name
	model_name=$(echo $model | sed -e 's/ /_/g')

	# pass env var to bsub
	bsub -env "all, model=$model_name" < run_deepbind_atac.lsf
done


