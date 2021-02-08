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

Nagano_meta = pd.read_csv("Nagano.meta.csv")
# ingest reference, get the working function
mm10_repli_score = funcs.repli_score("mm10_repli_chip.wig")

all_repli_score = funcs.collectcollect(Nagano_meta.query('passed_qc == 1'), "cell_nm", "pairs_0", mm10_repli_score, 128)

with open("out/all_repli_scores.pkl",'wb') as f:
    pickle.dump(all_repli_score,f)