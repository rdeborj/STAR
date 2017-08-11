from star.align import get_default_config, flatten_config_arrays, align_rnaseq

test_yaml = "test.yaml"
star_config = star.align.get_default_config(file=test_yaml)
star_config = flatten_config_arrays(config=star_config)
print(align_rnaseq(config=star_config))
