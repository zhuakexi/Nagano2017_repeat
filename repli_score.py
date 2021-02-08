import pandas as pd
import funcs
import pickle
from importlib import reload
reload(funcs)

# INPUT:
# threads: 128
# Nagano.meta.csv
# mm10_repli_chip.wig

# OUTPUT:
# repli_score.pkl

import ray
ray.init(address='auto',_redis_password='5241590000000000')
# read in datas sequencially
Nagano_meta = pd.read_csv("Nagano.meta.csv")
pairs_dict = {name:funcs.parse_pairs(other["pairs_0"]) for name, other in Nagano_meta.iloc[0:2,].set_index("cell_nm").iterrows()}
repli_chip = pd.read_table("mm10_repli_chip.wig", header=None, comment="#")
repli_chip.columns = "chr pos1 pos2 rep_value".split()
# tunning functions
f_ray_cell_repli_score = ray.remote(funcs.ray_cell_repli_score)
# do calc in cluster
res = [f_ray_cell_repli_score.remote(pairs_dict[key], repli_chip) for key in pairs_dict]
result = ray.get(res)
ray.shutdown()
with open("out/test.pkl",'wb') as f:
    pickle.dump(result,f)