from dataclasses import dataclass
from collections import OrderedDict
from typing import List
from distinctipy import get_colors, get_hex
from random import sample
import polars as pl
import polars.selectors as cs
@dataclass
class MetaInfo:
    meta_type: str
    vals: List[str]
    colors: List[str]
    name: str

@dataclass
class DimredInfo:
    name: str
    labels_2d: List[str]
    labels_3d: List[str]
    comps_2d: List[int]
    comps_3d: List[int]

@dataclass
class UIDefaults:
    default_subset_meta: str
    default_meta1: str
    default_meta2: str
    default_meta3: str
    default_dimred: str
    default_gene1: str
    default_gene2: str
    default_gene3: str
    default_genes_list: List[str]


class ScLiveDash:
    def __init__(self, adata, 
                 meta_ids=None,
                 meta_infos=None,
                 metas_schema=None,
                 dimred_ids=None,
                 default_dimred=None, 
                 default_subset_meta=None,
                 default_meta1=None, 
                 default_meta2=None, 
                 default_meta3=None, 
                 default_gene1=None, 
                 default_gene2=None, 
                 default_gene3=None, 
                 meta_colors=None,
                 meta_vals=None,
                 meta_names=None,
                 dimred_infos=None,
                 dimred_labels_2d = None,
                 dimred_labels_3d = None,
                 dimred_comps_2d = None,
                 dimred_comps_3d = None,
                 dimred_names = None,
                 default_genes_list=None):
        
        self.adata = adata
        self.meta_data = pl.from_pandas(adata.obs.reset_index()).rename({"index":"barcode"})
        self.meta_infos = OrderedDict()

        if meta_infos is None:
            meta_infos = OrderedDict()
        if metas_schema is None:
            metas_schema = OrderedDict()
        if meta_colors is None:
            meta_colors = OrderedDict()
        if meta_vals is None:
            meta_vals = OrderedDict()
        if meta_names is None:
            meta_names = OrderedDict()
        
        if meta_ids is None:
            meta_ids = self.meta_data.drop("barcode").columns
        
        for meta in meta_ids:
            if meta in meta_infos.keys():
                self.meta_infos[meta] = meta_infos[meta]
            else:
                self.add_meta(meta_id=meta,
                              meta_name=meta_names.get(meta, None),
                              meta_type=metas_schema.get(meta, None),
                              meta_vals=meta_vals.get(meta, None),
                              meta_colors=meta_colors.get(meta, None))
        
        if dimred_ids is None:
            self.dimred_infos = OrderedDict.fromkeys(self.adata.obsm.keys())
        else:
            self.dimred_infos = OrderedDict.fromkeys(dimred_ids)

        for dr in self.dimred_infos.keys():
            if dimred_infos and dr in dimred_infos.keys():
                self.dimred_infos[dr] = meta_infos[dr]
            else:
                if dimred_names and dr in dimred_names.keys():
                    name = dimred_names[dr]
                else:
                    name = dr.upper()
                
                if dimred_comps_2d and dr in dimred_comps_2d.keys():
                    comps_2d = dimred_comps_2d[dr]
                else:
                    comps_2d = [0,1]
                
                if dimred_comps_3d and dr in dimred_comps_3d.keys():
                    comps_3d = dimred_comps_3d[dr]
                else:
                    comps_3d = [0,1,2]
                
                if dimred_labels_2d and dr in dimred_labels_2d.keys():
                    labels_2d = dimred_labels_2d[dr]
                else:
                    labels_2d = [f"{name}_{i}" for i in range(2)]
                
                if dimred_labels_3d and dr in dimred_labels_3d.keys():
                    labels_3d = dimred_labels_3d[dr]
                else:
                    labels_3d = [f"{name}_{i}" for i in range(3)]
            self.dimred_infos[dr] = DimredInfo(name=name, labels_2d=labels_2d, labels_3d=labels_3d, comps_2d=comps_2d, comps_3d=comps_3d)
        if default_meta1 is None:
            default_meta1 = sample(list(self.meta_infos.keys()), 1)[0]
        if default_meta2 is None:
            default_meta2 = sample([k for k,v in self.meta_infos.items() if v.meta_type == "cat"], 1)[0]
        if default_meta3 is None:
            default_meta3 = sample([k for k,v in self.meta_infos.items() if v.meta_type == "cat" and k != default_meta2], 1)[0]
        if default_gene1 is None:
            default_gene1 = sample(list(self.adata.var_names), 1)[0]
        if default_gene2 is None:
            default_gene2 = sample(list(self.adata.var_names), 1)[0]
        if default_gene3 is None:
            default_gene3 = sample(list(self.adata.var_names), 1)[0]
        if default_subset_meta is None:
            default_subset_meta = sample([k for k,v in self.meta_infos.items() if v.meta_type == "cat"], 1)[0]
        if default_dimred is None:
            default_dimred = sample(list(self.dimred_infos.keys()), 1)[0]
        if default_genes_list is None:
            default_genes_list = sample(list(self.adata.var_names), 10)
        
        self.ui_defaults = UIDefaults(default_subset_meta=default_subset_meta, 
                                        default_meta1=default_meta1, 
                                        default_meta2=default_meta2, 
                                        default_meta3=default_meta3, 
                                        default_dimred=default_dimred, 
                                        default_gene1=default_gene1, 
                                        default_gene2=default_gene2, 
                                        default_gene3=default_gene3,
                                        default_genes_list=default_genes_list)
    def add_meta(self, meta_id,
                 meta_name=None,
                 meta_type=None,
                 meta_vals=None,
                 meta_colors=None):
        if meta_id in self.meta_infos.keys():
            raise ValueError(f"Meta {meta_id} already exists")
        elif meta_id not in self.meta_data.columns:
            raise ValueError(f"Meta {meta_id} not in adata.obs")
        else:
            if meta_name is None:
                meta_name = meta_id.upper()
            if meta_type is None:
                meta_type = "cat" if meta_id in self.meta_data.select(cs.string() | cs.categorical()).columns else "num"
            if meta_type == "cat":
                vals = self.meta_data[meta_id].unique().to_list()
                if meta_vals is None:
                    meta_vals = vals
                else:
                    if len(meta_vals) != len(vals):
                        raise ValueError(f"Meta values for {meta_id} are wrong length")
                if meta_colors is None:
                    meta_colors = [get_hex(c) for c in get_colors(len(vals))]
                else:
                    if len(meta_colors) != len(vals):
                        raise ValueError(f"Meta colors for {meta_id} are wrong length")
        self.meta_infos[meta_id] = MetaInfo(meta_type=meta_type, vals=meta_vals, colors=meta_colors, name=meta_name)
    def add_dimred(self, dimred_id,
                   dimred_name=None,
                   labels_2d=None,
                   labels_3d=None,
                   comps_2d=None,
                   comps_3d=None):
        if dimred_id in self.dimred_infos.keys():
            raise ValueError(f"Dimred {dimred_id} already exists")
        elif dimred_id not in self.adata.obsm.keys():
            raise ValueError(f"Dimred {dimred_id} not in adata.obsm")
        else:
            if dimred_name is None:
                dimred_name = dimred_id.upper()
            if labels_2d is None:
                labels_2d = [f"{dimred_name}_{i}" for i in range(2)]
            if labels_3d is None:
                labels_3d = [f"{dimred_name}_{i}" for i in range(3)]
            if comps_2d is None:
                comps_2d = [0,1]
            if comps_3d is None:
                comps_3d = [0,1,2]
        self.dimred_infos[dimred_id] = DimredInfo(name=dimred_name, labels_2d=labels_2d, labels_3d=labels_3d, comps_2d=comps_2d, comps_3d=comps_3d) 
    
    def remove_meta(self, meta):
        del self.sclive_dash_config["meta_info"][meta]

    def remove_metas(self, metas):
        for meta in metas:
            del self.sclive_dash_config["meta_info"][meta]

    def update_meta(self, **kwargs):
        for key, value in kwargs.items():
            if key in self.sclive_dash_config["meta_info"].keys():
                self.sclive_dash_config["meta_info"][key] = value
            else:
                raise ValueError(f"Meta {key} does not exist")
            
    def order_and_or_subset_metas(self, metas):
        self.sclive_dash_config["meta"] = OrderedDict([(meta, self.sclive_dash_config["meta_info"][meta]) for meta in metas])
