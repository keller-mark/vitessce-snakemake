__author__ = "Mark Keller"
__copyright__ = "Copyright 2022, Mark Keller"
__email__ = "mark_keller@hms.harvard.edu"
__license__ = "MIT"


import os
import json
from vitessce import VitessceConfig, MultiImageWrapper, OmeTiffWrapper, Component as cm

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

name = snakemake.params.get("name", "")
description = snakemake.params.get("description", "")
out_dir = snakemake.params.get("out_dir", "")

if hasattr(snakemake.input, "images"):
  images = snakemake.input.images
else:
  images = []

if hasattr(snakemake.input, "segmentations"):
  segmentations = snakemake.input.segmentations
else:
  segmentations = []

vc = VitessceConfig(name=name, description=description)
ds = vc.add_dataset(name=name).add_object(
  MultiImageWrapper(
    image_wrappers=([
      OmeTiffWrapper(img_path=image, name=f"Image ({i})")
      for i, image in enumerate(images)
    ] + [
      OmeTiffWrapper(img_path=image, name=f"Cell Segmentations ({i})", is_bitmask=True)
      for i, image in enumerate(segmentations)
    ]),
    use_physical_size_scaling=True,
  )
)
spatial = vc.add_view(cm.SPATIAL, dataset=ds)
lc = vc.add_view(cm.LAYER_CONTROLLER, dataset=ds)

vc.export(to="files", base_url="http://localhost:8000", out_dir=out_dir)

with open(snakemake.output[0], "w") as f:
  config_dict = vc.to_dict(base_url="http://localhost:8000")
  json.dump(config_dict, f)