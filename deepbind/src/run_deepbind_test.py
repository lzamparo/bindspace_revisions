import kipoi
import os
import sys

test = "D00290.003 # ALX3"

def process_id(line):
	line = line.strip()
	parts = line.split("#")
	code, tf_name = parts
	code = code.strip()
	tf_name = tf_name.lstrip()
	return code + "_SELEX_" + tf_name

prefix = "DeepBind/Homo_sapiens/TF/"
model_id = prefix + process_id(test)
print("Trying to get model: ", model_id)
try:
	model = kipoi.get_model(model_id)
except:
	print("Failed")

print("success hopefully")

