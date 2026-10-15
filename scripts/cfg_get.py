"""Print one value from the YAML config:  python3 cfg_get.py config.yaml paths.results_dir"""
import sys
import yaml

with open(sys.argv[1]) as fh:
    value = yaml.safe_load(fh)
for key in sys.argv[2].split("."):
    value = value[key]
print(value)
