import pandas as pd
import os
import gzip
from concurrent import futures
from functools import partial
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 
def window_count(distances:pd.DataFrame, win_num)->pd.Series:
    # count distribution of distance array
    windows = pd.Series([1000*2**(0.125*i) for i in range(0, win_num)]) # Peter's window
    window_count = []
    for index, point in enumerate(windows):
        if index == 0:
            count = len(distances[distances < point])
        elif index == win_num - 1:
            count = len(distances[distances >= point])
        else:
            count = len(distances[(distances >= point) & (distances < windows[index + 1])])
        window_count.append(count)
    window_count = pd.Series(window_count)
    window_count.index = range(1,win_num+1)
    # normalized by all intra contacts
    return window_count/len(distances)

def contact_describe(cell_name:str) -> pd.Series:
    # get cell's basic statistics, defined in Nagano2017
    contacts = parse_pairs(cell_name)
    intra = contacts.query(' chr1 == chr2 ')
    distances = abs(intra["pos1"] - intra["pos2"])
    
    all_ = len(distances[23_000 < distances])
    short = len(distances[(23_000 < distances) & (distances < 2_000_000)])
    mitotic = len(distances[(2_000_000 < distances) & (distances < 12_000_000)])
    farAvg = distances[(4_500_000 < distances) & (distances < 225_000_000)]
    
    mitotic_r = mitotic/all_
    short_r = short/all_
    
    # assign to different stages on Peter's cirtera
    if mitotic_r >= 0.3 and short_r <= 0.5:
        group = "Post-M"
    elif short_r > 0.5 and short_r + 1.8*mitotic_r > 1.0:
        group = "Pre-M"
    elif short_r <= 0.63:
        group = "G1"
    elif 0.63 < short_r <= 0.785:
        group = "early-S"
    elif short_r > 0.785:
        group = "late-S/G2"
    else:
        group = "blank"
    
    return pd.Series({"short%":short_r, "mitotic%":mitotic_r, "farAvg":farAvg.mean(),"group":group })
def collect(meta:pd.DataFrame, name_col:str, file_col:str, func:"callable", threads)->pd.DataFrame:
    # get cons for cells in dir
    file_names = meta[file_col]
    with futures.ProcessPoolExecutor(threads) as pool:
        res = pool.map(func, file_names)
    result = pd.concat(res, axis=1)
    result.columns = meta[name_col]
    return result
def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #read table format data
    pairs = pd.read_table(filename, header=None, comment="#",low_memory=False)
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # get real sample name
    #assign column names
    pairs.columns = name_array[0:pairs.shape[1]]
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def dis_counts(cell_name:str):
    # work for 11 column table only
    # get cell's intra contact's distribution in Peter's window
    # using customized .pairs parser
    contacts = parse_pairs(cell_name)

    # get contact distance array
    intra = contacts.query("chr1 == chr2")
    distances = abs(intra["pos1"] - intra["pos2"])
    # count according to Peter's window
    counts = window_count(distances, 150)
    counts.name = cell_name
    #return counts
    return counts.reindex(range(38,151)) # only show 38-150
# calculate single cell replicate score 
def repli_score(repli_file:str):
    # need wig or bed file
    return partial(cell_repli_score, repli_file=repli_file)
def get_mean_repli_value(line:pd.Series, full_ref:pd.DataFrame)->int:
    # line: one leg
    # table: full repli_ref
    chr_name, pos_range = line.name[0], line.name[1]
    by_chr_ref = full_ref.loc[chr_name]
    hit = by_chr_ref[by_chr_ref.index.overlaps(pos_range)] # line name: (chr, interval), use interval
    if len(hit) == 0:
        return pd.NA
    else:
        return hit["rep_value"].mean()
def chunck_repli_values(chunck:pd.DataFrame, full_ref:pd.DataFrame)->pd.DataFrame:
    # transform pairs chunck to repli-value chunck
    return chunck.apply(get_mean_repli_value, full_ref=full_ref, axis=1)
def cell_repli_score(pairs_file:str, repli_file:str)->dict:
    # read reference file and target pairs file
    repli_chip = pd.read_table(repli_file, header=None, comment="#")
    repli_chip.columns = "chr pos1 pos2 rep_value".split()
    pairs = parse_pairs(pairs_file)

    # get legs from contacts
    leg1, leg2 = pairs[["chr1","pos1"]], pairs[["chr2","pos2"]]
    leg1.columns = leg2.columns = "chr pos".split()
    legs = leg1.append(leg2)
    # expand 10_000 each way
    legs_lower, legs_upper = legs["pos"] - 10_000, legs["pos"] + 10_000
    legs_lower[legs_lower < 0] = 0
    legs_lower.name, legs_upper.name = "lower","upper"
    legs = pd.concat([legs["chr"],legs_lower, legs_upper],axis=1)

    # using pd.IntervalIndex to do searching
    ## build index for reference
    chip_range_index = pd.IntervalIndex.from_arrays(repli_chip["pos1"],repli_chip["pos2"])
    chip_range_index.name = "range_index"
    repli_chip_indexed = repli_chip.set_index([repli_chip["chr"], chip_range_index])
    ## build index for target
    legs_range_index = pd.IntervalIndex.from_arrays(legs["lower"],legs["upper"])
    legs_range_index.name = "pos_range"
    legs_indexed = legs.set_index([legs["chr"],legs_range_index])

    # do calc(time consuming)
    results = {}
    for key, per_chr_legs in legs_indexed.groupby(level=0):
        results[key] = chunck_repli_values(per_chr_legs, repli_chip_indexed) # return series

    # tidy
    result_frame = pd.DataFrame()
    for key in results:
        result_frame = result_frame.append(results[key].reset_index()[["chr",0]], ignore_index=True)
    result_frame.columns = "chr value".split()

    raw_fends = len(result_frame)
    annote_fends = len(result_frame.dropna())
    early_replicate_fends = len(result_frame.query('value >0 '))  

    return {"repli_score":early_replicate_fends/annote_fends, "annote_ratio":annote_fends/raw_fends}