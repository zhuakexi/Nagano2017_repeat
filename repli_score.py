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
Nagano_meta = pd.read_csv("Nagano.meta.csv")
# remote version of working function
f_cell_repli_score = ray.remote(funcs.cell_repli_score)
res = [f_cell_repli_score.remote(file_name, "mm10_repli_chip.wig") for file_name in Nagano_meta["pairs_0"][0:2]]
result = ray.get(res)
ray.shutdown()

with open("out/test.pkl",'wb') as f:
    pickle.dump(result,f)