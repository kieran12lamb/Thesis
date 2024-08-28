import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
import re
import json
from sklearn.utils import resample

#The ClusterViz class is used to generate a number of plots for mutational signature analysis.
class SignaturePreprocessor:
    #Constructor for the cluster visualisation class
    def __init__(self, data, metadata):
        self.mutational_burden = data.sum(axis=1)
        self.metadata = metadata.dropna(axis=1, how='all')

        self.data = data 
        mask = data.columns.str[2] !="X"
        self.data = self.data[data.columns[mask]]
        self.metadata = metadata

         # drop rows with no lineage
        nan_lineages = self.metadata.lineage[self.metadata.lineage.isna()].index
        self.metadata = self.metadata.drop(nan_lineages)
        self.data = self.data.drop(nan_lineages)

        if "epi_week" in self.metadata.columns:
            nan_epi_weeks= self.metadata.epi_week[self.metadata.epi_week.isna()].index
            self.metadata = self.metadata.drop(nan_epi_weeks)
            self.data = self.data.drop(nan_epi_weeks)
            self.metadata.epi_week = self.metadata.epi_week.astype(int)
        self.metadata = self.metadata.replace(["NA","NAN","NaN"],"")
        
    def normalisations(self,percentages=False,reduce_resample=False,clustering=[0,None],pseudosample_on=None, resample=None,reduce=None,min_max_norm = False):
        if resample is not None:    
            self.data = self.resample(self.data, resample)
        elif reduce is not None:
            self.data = self.reduction(self.data)
        elif reduce_resample == True:
            self.data = self.reduce_resample(samples=None,pseudosample_on= pseudosample_on)
        if clustering[0] != 0:
            self.cluster(num_clusters=clustering[0], method=clustering[1])
            
        if min_max_norm == True:
            self.data=(self.data-self.data.min())/(self.data.max()-self.data.min())
        
        if pseudosample_on is not None:
            if type(pseudosample_on) is str:
                self.data = self.data.groupby(by=self.metadata[pseudosample_on]).sum()
                self.metadata = self.make_pseudosample_metadata(pseudosample_on= pseudosample_on)
            else:
                self.pseudo_size = self.data.groupby(by=[self.metadata[column] for column in pseudosample_on]).size().reset_index()
                self.pseudo_size.index = self.data.groupby(by=[self.metadata[column] for column in pseudosample_on]).size().reset_index().astype(str)[pseudosample_on].apply('_'.join, axis=1)
#                 self.pseudo_size.index =pd.Series(self.pseudo_size[self.pseudo_size.columns[0]],dtype=str) + "_" + pd.Series(self.pseudo_size[self.pseudo_size.columns[1]],dtype=str)
                self.pseudo_size = self.pseudo_size.drop(pseudosample_on,axis=1)
    
                data_index =  self.data.groupby(by=[self.metadata[column] for column in pseudosample_on]).size().reset_index().astype(str)[pseudosample_on].apply('_'.join, axis=1)
                self.data = self.data.groupby(by=[self.metadata[column] for column in pseudosample_on]).sum().reset_index()
                self.data.index = data_index
                self.metadata = self.make_pseudosample_metadata(pseudosample_on= pseudosample_on)
        #Convert data into percentages
        if percentages == True:
            self.sum_mutations = self.data.sum(axis=1)
            # data = self.data.replace(0, np.nan)
            data = self.data
            data = data.div(data.sum(axis=1),axis=0)
            self.data = data.fillna(0)
        #Keep data as raw counts
        else:
            self.sum_mutations = self.data.sum(axis=1)
        
        self.data = self.data.sort_index()
        self.metadata = self.metadata.sort_index()

    #Cluster the data using either kmeans or gmm clustering algorithm
    def cluster(self,num_clusters,method = "kmeans", seed=None):
        data = self.data
        if method == "kmeans":
            kmeans = KMeans(n_clusters=num_clusters,random_state=seed).fit(data)
            clusters = kmeans.predict(data)
        elif method == "gmm":
            gmm = GaussianMixture(n_components=num_clusters,random_state=seed).fit(data)
            clusters = gmm.predict(data)   
        self.metadata['Cluster'] = clusters

    def reduce_resample(self,samples=None,pseudosample_on=None):
        V = self.data
        V = V.sample(n=len(V), replace=True)
        V = self.reduction(V)
        V = V.fillna(0)
        return V
            
    def reduction(self,V):
        summed_V = V.sum(axis=0,numeric_only=True).values/V.sum(axis=0,numeric_only=True).sum(axis=0)
        summed_V = pd.Series(index=V.sum(axis=0,numeric_only=True).index,data=summed_V)
        reduced_V = pd.DataFrame(V[summed_V.index[summed_V > 0.005]])
        full = pd.DataFrame(columns =V.columns).drop(reduced_V.columns, axis=1,)
        full = pd.merge(reduced_V,full,right_index=True,left_index=True,how="left" )
        full = full.fillna(0)
#         full = V
        return full

    def resample(self,V,samples):
        if samples != None:
            V = resample(V, replace=True, n_samples=samples)
        else:
            V = resample(V, replace=True, n_samples=len(V))
        shape = V.shape
        V = pd.merge(V,self.metadata, left_index=True,right_index=True,how="left")
        self.metadata = V[V.columns[shape[1]:]]
        V = V[V.columns[:shape[1]]]
        return V


    def make_pseudosample_metadata(self,pseudosample_on):
        if isinstance(pseudosample_on,str) == True:
            pseudo_groups = self.metadata[pseudosample_on].unique()
            
            pseudo_metadata = []
            for group in pseudo_groups:
                # pseudo_group = self.metadata.drop(['sample_date','lineage_support','pillar_2',],axis=1)[self.metadata[pseudosample_on] == group]
                pseudo_group = self.metadata.drop(["GISAID_ID","Host","Non-Shortcut-Lineage","Length","sample_date","FullCountry","Error"],axis=1)[self.metadata[pseudosample_on] == group]
                # pseudo_group['Cluster'] = pseudo_group['Cluster'].apply(str)
                pseudo_group['epi_week'] = pseudo_group['epi_week'].apply(str)
                num_sequences = len(pseudo_group)
                cat_metadata = pd.get_dummies(pseudo_group)
                cat_metadata = cat_metadata.sum(axis=0)

                cat_metadata = cat_metadata.fillna(0)
                cat_metadata['num_sequences'] = num_sequences
                pseudo_metadata.append(cat_metadata)

            pseudo_metadata = pd.concat(pseudo_metadata,axis=1).T
            pseudo_metadata.index = pseudo_groups
            pseudo_metadata.fillna(0).sort_index()
            # pseudo_metadata = pseudo_metadata.drop(pseudosample_on)
        else:
            pseudo_metadata = self.data[pseudosample_on]
            self.data = self.data.drop(pseudosample_on,axis=1)
            if "epi_week" in pseudo_metadata.columns:
                pseudo_metadata['epi_week'] = pd.Series(pseudo_metadata["epi_week"],dtype=int)
            else:
                pseudo_metadata['sample_date'] = pd.Series(pseudo_metadata["sample_date"],dtype=int)
            pseudo_metadata["Size"] = self.pseudo_size[0]
            pseudo_metadata
        return pseudo_metadata

