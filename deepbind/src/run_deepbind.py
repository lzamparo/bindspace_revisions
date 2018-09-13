import kipoi
import os
import sys
import time
from dataloader import TruncSeqDataset as DL


def process_id(line):
	line = line.strip()
	parts = line.split("_#_")
	code, tf_name = parts
	code = code.strip()
	tf_name = tf_name.strip()
	return code + "_SELEX_" + tf_name

def generate_output_file(line):
	line = line.strip()
	parts = line.split("_#_")
	_, tf_name = parts
	return tf_name + "_outputs.txt"

def determine_human_vs_mouse(line):
	line = line.strip()
	parts = line.split("_#_")
	_, tf_name = parts
	if tf_name.isupper():
		return "DeepBind/Homo_sapiens/TF/"
	else:
		return "DeepBind/Mus_musculus/TF/"

# Parse argv to get the model specification
tf_line = sys.argv[1]

# Bookeeping for input, output resolution
model_prefix = determine_human_vs_mouse(tf_line)
data_prefix = "~/projects/revisions/bindspace_revisions/deepbind/data"
output_prefix = "~/projects/revisions/bindspace_revisions/deepbind/outputs"

k562_atlas_bed = os.path.join(os.path.expanduser(data_prefix), "K562_atac.bed")
k562_atlas_fasta = os.path.join(os.path.expanduser(data_prefix), "fixed_atac.fasta")

model_id = model_prefix + process_id(tf_line)
outfile_name = generate_output_file(tf_line)
output_file = os.path.join(os.path.expanduser(output_prefix), outfile_name)

# Fetch the model with kipoi
print("Trying to get model: ", model_id, flush=True)
try:
	model = kipoi.get_model(model_id)
except:
	print("Failed to process model_id")

print("Loaded the model", flush=True)

# setup the example dataloader kwargs
dl_kwargs = {'intervals_file': k562_atlas_bed, 'fasta_file': k562_atlas_fasta}

# Option 1: use data loader instantiation of an iterator directly
dl = DL(**dl_kwargs)

t0 = time.time()
with open(output_file, 'w') as outfile:
	batch_iter = dl.batch_iter(batch_size=128)
	for batch in batch_iter:
		preds = model.predict_on_batch(batch['inputs'])
		for p in preds:
			print(p, file=outfile)
t1 = time.time()

print("Took ", t1 - t0, "s to complete", flush=True)

# Option 2: use pipeline.predict passing the dl_kwargs dict directly 
# for prediction on all sites in the atlas in stages  

#t0 = time.time()
#try:
#	preds_on_atlas = model.pipeline.predict(dl_kwargs, batch_size=128)
#except AssertionError as ae:
#	print("Assertion error of some kind:", ae)
#t1 = time.time()
#print("Took ", t1 - t0, "s to complete", flush=True)

#with open(outfile_name, 'w') as outfile:
#	for p in preds_on_atlas:
#		print(p, file=outfile)

