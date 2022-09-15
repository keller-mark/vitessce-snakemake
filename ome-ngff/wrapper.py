__author__ = "Mark Keller"
__copyright__ = "Copyright 2022, Mark Keller"
__email__ = "mark_keller@hms.harvard.edu"
__license__ = "MIT"


import os
from os.path import join
import json
import zarr
from vitessce import VitessceConfig, AbstractWrapper, Component as cm

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

name = snakemake.params.get("name", "")
description = snakemake.params.get("description", "")
out_dir = snakemake.params.get("out_dir", "")
port = snakemake.params.get("port", "")

if hasattr(snakemake.input, "images"):
  images = snakemake.input.images
else:
  images = []

class OmeNgffWrapper(AbstractWrapper):
  def __init__(self, store, **kwargs):
    super().__init__(**kwargs)
    self._store = store
    self.zarr_folder = 'image.ome.zarr'

  def convert_and_save(self, dataset_uid, obj_i):
    super().convert_and_save(dataset_uid, obj_i)
    zarr_filepath = self.get_zarr_path(dataset_uid, obj_i)
    dest = zarr.DirectoryStore(zarr_filepath)
    zarr.copy_store(self._store, dest)

    self.file_def_creators += [
      self.make_raster_file_def_creator(dataset_uid, obj_i)
    ]
    self.routes += self.get_out_dir_route(dataset_uid, obj_i)

  def get_zarr_path(self, dataset_uid, obj_i):
    out_dir = self._get_out_dir(dataset_uid, obj_i)
    zarr_filepath = join(out_dir, self.zarr_folder)
    return zarr_filepath

  def get_zarr_url(self, base_url="", dataset_uid="", obj_i=""):
    return self._get_url(base_url, dataset_uid, obj_i, self.zarr_folder)
  
  def make_raster_file_def_creator(self, dataset_uid, obj_i):
    def get_raster(base_url):
      obj_file_def = {
        "type": "raster",
        "fileType": "raster.ome-zarr",
        "url": self.get_zarr_url(base_url, dataset_uid, obj_i)
      }
      return obj_file_def
    return get_raster

vc = VitessceConfig(name=name, description=description)
ds = vc.add_dataset(name=name)
for image in images:
  store = zarr.DirectoryStore(image)
  ds.add_object(OmeNgffWrapper(store=store))
spatial = vc.add_view(cm.SPATIAL, dataset=ds)
lc = vc.add_view(cm.LAYER_CONTROLLER, dataset=ds)

vc.layout(spatial | lc)

vc.export(to="files", base_url=f"http://localhost:{port}", out_dir=out_dir)

with open(snakemake.output[0], "w") as f:
  config_dict = vc.to_dict(base_url=f"http://localhost:{port}")
  json.dump(config_dict, f)

print("To open Vitessce, serve the processed files and then open in a web browser:")
print(f"http-server --cors='*' --port {port} ./data/")
print(f"open \"http://dev.vitessce.io/?url=http://localhost:{port}/vitessce.config.json\"")