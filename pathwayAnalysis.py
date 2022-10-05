import pandas as pd
import sys
sys.path.append('../GeneralMethods/')
from statsmodels.stats.multitest import fdrcorrection
import scipy.stats as ss
import os
import re


def RobustZScore(arr, lb=10**(-5), ub=30):
    
    ind, = np.nonzero(np.abs(arr)<lb)
    if len(ind)>0:
        arr[ind] = np.zeros(len(arr[ind]))
    
    STD = np.std(arr) if np.std(arr)!=0 else 1
    res = (arr-np.mean(arr))/STD
    ind2, = np.nonzero(np.abs(res)>ub)
    
    if len(ind2)>0:
        sign = arr[ind2]/np.abs(arr[ind2])
        ind3, = np.nonzero(np.abs(res)<=ub)
        tmp = res[ind3]
        res[ind2] = sign*np.max(np.abs(tmp))
    
    if len(arr)==len(res):
        return res


def read_met_model(path):
    import scipy.io
    mat = scipy.io.loadmat(path)
    mat = mat[[k for k in mat.keys()][-1]]
    keys = mat.dtype.fields.keys()
    arrs = mat.item()
    model_dict = {k:arr for k, arr in zip(keys, arrs)}
    return model_dict


def hypergeom_test(arr1, arr2, label=1):
    """
    arr1: array with only 1, 0, -1
    arr2: array with only 1, 0, -1
    """
    from scipy.stats import hypergeom
    M = len(arr1)
    N = sum(arr1==label)
    n = sum(arr2==label)
    k = sum((arr1==label)*(arr2==label))
    # more overlaps higher pvalues
    hypergeom_p = hypergeom.sf(k-1, M, n, N)
    a, b, c, d = k, n-k, N-k, M-(n+N)+k
    table = np.array([[a, b], [c, d]])
    start, end = hypergeom.support(M, n, N)
    # print(sum(hypergeom.pmf(np.arange(k, end+1), M, n, N)))
    from scipy.stats import fisher_exact
    oddsr, fisher_p = fisher_exact(table, alternative='greater')
    # print(k-1, M, n, N, fisher_p, hypergeom_p)
    return hypergeom_p, fisher_p#, M, N, n, k

def ranksumtest(arr1, arr2):
    from scipy.stats import ranksums
    w, p = ranksums(arr1, arr2, alternative='less')
    return w, p




def pathwayAnalysis(model_path='./model_MGSA.mat', subSysOfInterest=[['all']]):

    import re
    import numpy as np
    
    # load metabolic model
    recon1 = read_met_model('/home/daweilin/StemCell/model_MGSA.mat')
    recon1g = [str(arr[0])[2:-2] for arr in recon1['genes']]
    recon1n = [str(arr) for arr in recon1['ngenes']]
    
    # subsystems
    subSysIter = np.unique(recon1['subSystems']) if subSysOfInterest==[['all']] else subSysOfInterest
    
    # group reactions to belonging pathways/subsystems
    subsys_genes = {k[0]:[] for k in subSysIter}

    # iterate thru all interested subsys
    for subsys in subSysIter:
        subsys = subsys[0]
        belonging_genes = []
        for arr in recon1['grRules'][recon1['subSystems']==subsys]:
            if len(arr)>0:
                s = arr[0]
                for ele in re.findall('\d+\.\d+', s):
                    belonging_genes.append(ele)
        belonging_genes = [ele.split('.')[0] for ele in belonging_genes]
        
        # convert gene IDs to symbols
        import mygene
        mg = mygene.MyGeneInfo()
        query = mg.getgenes(np.unique(belonging_genes),
                    fields=['symbol', 'entrezgene'])
        
        # gather genes sumbol
        gene_sym_arr = []
        for sym in query:
            try:
                gene_sym_arr.append(sym['symbol'])
            except:
                print(sym, 'not found')
        
        subsys_genes[subsys] = gene_sym_arr

    return subsys_genes


def enrichmentAnalysis(input_df_path, colname_to_ID, subsys_dict, method='ranksum'):
    
    # method seleciton
    testMethod = ranksumtest if method=='ranksum' else hypergeom_test

    # import data
    input_df = pd.read_csv(input_df_path)
    
    subsys_pv = {k:[] for k in subsys_dict.keys()}
    # iterate thru the columns (different cells)
    for col in input_df.columns:
        if col!=colname_to_ID:
            for k in subsys_dict.keys():
                sym_arr = subsys_dict[k]
                scores = RobustZScore(input_df[col].to_numpy())
                _, p = testMethod(
                        scores[input_df[colname_to_ID].isin(sym_arr)],
                        scores[input_df[colname_to_ID].isin(sym_arr)==0]
                        )
                subsys_pv[k].append(p)

    return subsys_pv



def pathwayAnalysisReaction(model_path='./model_MGSA.mat', subSysOfInterest=[['all']])

    import re
    import numpy as np
    
    # load metabolic model
    recon1 = read_met_model('/home/daweilin/StemCell/model_MGSA.mat')
    recon1g = [str(arr[0])[2:-2] for arr in recon1['genes']]
    recon1n = [str(arr) for arr in recon1['ngenes']]

    # subsystems
    subSysIter = np.unique(recon1['subSystems']) if subSysOfInterest==[['all']] else subSysOfInterest
    
    # group reactions to belonging pathways/subsystems
    subsys_rxns = {k[0]:[] for k in subSysIter}

    # iterate thru all interested subsys
    for subsys in subSysIter:
        subsys = subsys[0]
        belonging_rxns = []
        for ele in recon1['rxns'][recon1['subSystems']==subsys]:
            belonging_rxns.append(ele)
        
        subSys_rxns[subsys] = belonging_rxns

    return subSys_rxns


