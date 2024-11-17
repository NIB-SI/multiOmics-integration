# python modules
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from scipy.signal import medfilt
from scipy.signal import argrelextrema
from math import isclose
import glob

import os
import seaborn as sns
import matplotlib.pyplot as plt

def get_interative_t_values(F, D, d_min_t=None):
    
    if d_min_t is None:
        d_min_t = {"general":0}
    
    all_dfs = []   
    for f in F:
        print(f)
        df = pd.read_csv(f, sep="\t")    
        
        if df.shape[0] == 0:
            continue
        
        else:
            all_dfs.append(df)
        
    if len(all_dfs) == 0:
        return pd.DataFrame()
    
    df = pd.concat(all_dfs).groupby("threshold").max()# skipna=True)    
    df.sort_values("threshold", inplace=True)
    df.reset_index(inplace=True)
    df = df[df["threshold"] > d_min_t["general"]]
    df.reset_index(inplace=True, drop=True)

    D['aplestin'] = np.nan
    D['cc_inflection'] = np.nan
    D['density_min'] = np.nan
    D['elo_clustering'] = np.nan
    D['gupta_clustering'] = np.nan
    D['single_component'] = np.nan
    D['whole_graph'] = np.nan        

    if df.shape[0] <3:
        return df
    
    # single component
    # First point before more than one cc appears
    more_one_cc = df[df["connected-component-count"] > 1]
    if more_one_cc.shape[0] > 0:
        prev_arg = more_one_cc.iloc[0].name
        if prev_arg > 1:
            D['single_component'] = df["threshold"][prev_arg -1]
    
    # whole graph
    D['whole_graph'] = df[df["vertex-count"] ==  df["vertex-count"].max()]["threshold"].max()

    # giant cc down inflection
    diffs = np.diff(df["largest-cc-size"].ewm(span = 20).mean())
    drop_cutoff = diffs.std()
    diffs_b = (abs(diffs) > drop_cutoff) & (diffs < 0)
    if sum(diffs_b) > 0:
        D['cc_inflection'] = df["threshold"][np.argmax(diffs_b) + df.index[0]]
        
    # density minimum
    D['density_min'] = df[df["density"] ==  df["density"].min()]["threshold"].min()   
        
    # apelstin
    Nsv = df["edge-count"]/df["vertex-count"]
    dNsv_dt = np.gradient(Nsv, df["threshold"])
    df["dNsv_dt"] = dNsv_dt
    D['aplestin'] = df["threshold"][dNsv_dt > 0].min()

    if not np.all(np.isnan(df["clustering-coefficient"])):
        # gupta clustering coefficient
        # largest t with sharp increase in clustering-coefficient
        # estimate by first point where at least 
        # 3 differences are larger than 0.3 * stddev
        found_gupta = False
        diffs = np.diff(df["clustering-coefficient"])
        if (np.sum(diffs>0) > 3): # first check: is there consistent increase in clustering coeficeint?
            cutoff = np.nanstd(diffs[diffs>0])*0.3
            prev_prev_d = diffs[0]
            prev_d = diffs[1]
            i = 0
            for d in diffs[2:]:
                if np.all(np.array([prev_prev_d, prev_d, d]) > cutoff) and (i > 1):
                    found_gupta = True
                    break
                prev_prev_d = prev_d
                prev_d = d
                i += 1

        if found_gupta:
            D['gupta_clustering'] = df["threshold"][i]
        
        if not np.all(np.isnan(df["random-clustering-coefficient"])):
            # elo clustering coefficient
            # first local maximum
            C0_diffs =  medfilt(df["clustering-coefficient"] - df["random-clustering-coefficient"])
            df["elo-clustering-coefficient-diffs"] = C0_diffs
            D['elo_clustering'] = df["threshold"][argrelextrema(C0_diffs, np.greater_equal)[0][0] + df.index[0]].min()

    return df


def get_result_per_prefix(result_path, prefix):
    iterative_results = result_path.glob(prefix + ".[0-9]*.iterative.txt")
    D = {}
    df = get_interative_t_values(iterative_results, D)
    
    return D, df
