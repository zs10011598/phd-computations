import os
import psycopg2
import datetime
import psycopg2.extras
from sqlalchemy import create_engine
import pandas as pd
import numpy as np


def get_results_covariable(data_train, df_covocc, df_cov, data_pob_sample):
    """
        Description: 
    """
    N = data_pob_sample['pob'].sum()
    Nc = data_train['cases'].sum()
    PC = Nc/N
    Nc_ = N - Nc
    P_C = 1 - PC
    s0 = np.log(PC/P_C)
    alpha = 0.01
    
    dict_results = {
                    'name': [],
                    'id': [],
                    'variable': [], 
                    'Nx': [], 
                    'Ncx': [], 
                    'PCX': [], 
                    'PC': [], 
                    'Nc': [], 
                    'N': [], 
                    'epsilon': [],
                    'Nc_': [],
                    'Nc_x':[],
                    'P_C': [],
                    'P_CX':[],
                    's0':[],
                    'score':[]}

    for index, row in df_cov.iterrows():
        
        dict_results['name'].append(row['name'] + ' ' + row['interval'])
        dict_results['variable'].append(row['code'] + ' ' + row['interval'])
        dict_results['id'].append(row['id'])
        
        Nx = 0
        Ncx = 0
        
        df = df_covocc[df_covocc['covariable_id'] == row['id']]
        occurrence_covariable_list = df['gridid_mun'].tolist()
        df = None
        Ncx = data_train[data_train['gridid_mun'].isin(occurrence_covariable_list)]['cases'].sum()
        Nx = data_pob_sample[data_pob_sample['gridid_mun'].isin(occurrence_covariable_list)]['pob'].sum()
        
        dict_results['Nx'].append(Nx)
        dict_results['Ncx'].append(Ncx)
        
        if Nx == 0:
            PCX = 0
        else:
            PCX = Ncx/Nx
        
        dict_results['PCX'].append(PCX)
        dict_results['PC'].append(PC)
        dict_results['Nc'].append(Nc)
        dict_results['N'].append(N)
        
        if Nx == 0 or PC == 0:
            epsilon = 0
        else:
            epsilon = (Nx*(PCX - PC))/ ((Nx*PC*(1 - PC))**0.5)

        dict_results['epsilon'].append(epsilon)
        dict_results['Nc_'].append(Nc_)
        
        Nc_x = Nx - Ncx
        dict_results['Nc_x'].append(Nc_x)
        
        dict_results['P_C'].append(P_C)
        
        if Nx == 0:
            P_CX = 1
        else:
            P_CX = Nc_x/Nx
    
        dict_results['P_CX'].append(P_CX)
        
        dict_results['s0'].append(s0)
        
        score = np.log(((Ncx + 0.005)/(Nc + 0.01))/((Nc_x + 0.01)/(Nc_+0.005)))
        dict_results['score'].append(score)
        
    return pd.DataFrame(dict_results)


def get_results_cell(df_results_cov_training, data_train, df_covocc, data_pob_sample):
    """
        Description:
    """
    scores_municipalities_training = {}
    s0 = df_results_cov_training.iloc[0]['s0']

    for index, row in data_pob_sample.iterrows():
        scores_municipalities_training[row['gridid_mun']] =  {'score': s0, 'gridid_mun': row['gridid_mun'], 'pob': row['pob']}
    
    for index, row in df_results_cov_training.iterrows():
        municipalities = df_covocc[df_covocc['covariable_id'] == row['id']]['gridid_mun'].tolist()
    
        for municipality in municipalities:
            scores_municipalities_training[municipality]['score'] += row['score']
    
    df_cells = pd.DataFrame(scores_municipalities_training.values())
    df_cells = df_cells.merge(data_train, how='left', on='gridid_mun').sort_values('score', ascending=False)
    df_cells['cases'] = df_cells['cases'].fillna(0)
    df_cells = df_cells.reset_index(drop=True)

    N = data_pob_sample['pob'].sum()
    percentil_length = N / 20
    aux_training = 0
    cummulated_length = 0
    p_training = []
    percentil_training = []

    for d in range(20):
        lower_training = aux_training
        upper_training = aux_training
    
        while cummulated_length < (d+1)*percentil_length:    
            cummulated_length += df_cells.iloc[upper_training]['pob']
            upper_training += 1
    
        aux_training = upper_training
        cases_percentil_training = df_cells.iloc[lower_training:upper_training]['cases'].sum()
        pobtot_percentil_training = df_cells.iloc[lower_training:upper_training]['pob'].sum()
        p_training += [cases_percentil_training/pobtot_percentil_training for i in range(upper_training - lower_training)]
        percentil_training += [20 - d for i in range(upper_training - lower_training)]
    
    df_cells['p'] = pd.Series(p_training)
    df_cells['percentile'] = pd.Series(percentil_training)

    return df_cells