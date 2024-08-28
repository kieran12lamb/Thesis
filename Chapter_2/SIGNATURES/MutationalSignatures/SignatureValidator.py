import numpy as np
import pandas as pd
import re
import json
import MutationalSignatures.SignaturePreprocessor as SignaturePreprocessor
import MutationalSignatures.SignatureVisualiser as SignatureVisualiser
import MutationalSignatures.compiled_functions as cf

from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn.metrics.pairwise import cosine_similarity
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.metrics import pairwise_distances

class SignatureValidator:
    #Constructor for the signature validator class
    def __init__(self, data, metadata,percentages=False,clustering=[0,None],pseudosample_on=None,cluster_distance="cosine",sequence_filter = 0,traditional=False):
        # Sig Vis Variables
        self.data = data
        self.metadata = metadata
        
        self.sequence_filter = sequence_filter
    
        # Sig Pre-processing variables
        self.percentages=percentages
        self.clustering=clustering
        self.pseudosample_on=pseudosample_on
        self.cluster_distance=cluster_distance
        self.traditional=traditional
        
        self.sig_process = SignaturePreprocessor.SignaturePreprocessor(data, metadata)
        self.sig_process.normalisations(percentages=percentages,clustering=clustering,pseudosample_on=pseudosample_on,reduce=True)
        
        # Mutational Catalogue Replicates
        self.mutational_catalogues = {}
        self.mutational_exposures = {}
        self.error = {}
        
    def construct_replicates(self,MIN_N = 2, MAX_N = 15, iterations = 10):
        self.MIN_N = MIN_N
        self.MAX_N = MAX_N
        for i in range(MIN_N,MAX_N):
            catalogue_replicates = []
            exposure_replicates = []
            error_replicates = []
            for j in range(iterations):
                sig_process = SignaturePreprocessor.SignaturePreprocessor( 
                        data=self.data,              
                        metadata=self.metadata)
                print(f"N={i}, iteration={j}")
                sig_process.normalisations(percentages= self.percentages,
                                            reduce_resample= True,
                                            clustering=self.clustering,
                                            pseudosample_on=self.pseudosample_on
                                            )
                if self.pseudosample_on is not None:
                    sig_vis = SignatureVisualiser.SignatureVisualiser(sig_process.data[sig_process.metadata.Size>=self.sequence_filter],sig_process.metadata[sig_process.metadata.Size>=self.sequence_filter],i,self.pseudosample_on,traditional=self.traditional,)
                    catalogue_replicates.append(sig_vis.signatures.T)
                    new_exposure = pd.DataFrame(sig_vis.exposures)
                    new_exposure.columns = sig_process.data[sig_process.metadata.Size>=self.sequence_filter].index
                    exposure_replicates.append(new_exposure)
                    error_replicates.append(sig_vis.nmf.reconstruction_err_)
                else:
                    sig_vis = SignatureVisualiser.SignatureVisualiser(sig_process.data,sig_process.metadata,i,self.pseudosample_on,traditional=False,)
                    catalogue_replicates.append(sig_vis.signatures.T)
                    new_exposure = pd.DataFrame(sig_vis.exposures)
                    new_exposure.columns = sig_process.data.index
                    exposure_replicates.append(new_exposure)
                    error_replicates.append(sig_vis.nmf.reconstruction_err_)
            self.mutational_catalogues[str(i)] = catalogue_replicates
            self.mutational_exposures[str(i)] = exposure_replicates
            self.error[str(i)] = error_replicates
                
    def metrics(self):
        scores = {"Silhouette Score for N Signatures":[],"Reconstruction Error for N Clusters":[],"Silhouette Score Per Cluster":[]}
        for i in range(self.MIN_N,self.MAX_N):
            catalogue = pd.concat(self.mutational_catalogues[str(i)],axis=0).reset_index(drop=True)

            distances = pairwise_distances(catalogue,catalogue,metric=self.cluster_distance)
            clusters = KMeans(n_clusters=i).fit_predict(distances)
            scores["Silhouette Score for N Signatures"].append(metrics.silhouette_score(X=distances, labels=clusters, metric="precomputed"))
            
            per_sample_silhouette = metrics.silhouette_samples(distances, clusters, metric="precomputed")
            cluster_silhouette = {}
            for idx, label in enumerate(clusters):
                if str(label) in cluster_silhouette:
                    cluster_silhouette[str(label)].append(per_sample_silhouette[idx])
                else:
                    cluster_silhouette[str(label)] = [per_sample_silhouette[idx]]
            for label in cluster_silhouette.keys():
                cluster_silhouette[label] = [np.mean(cluster_silhouette[label])]
            cluster_silhouette = pd.DataFrame(cluster_silhouette).T.sort_index()
            cluster_silhouette.columns = [str(i)]
            scores["Silhouette Score Per Cluster"].append(cluster_silhouette)
            scores["Reconstruction Error for N Clusters"].append(np.mean(self.error[str(i)]))
        scores["Silhouette Score Per Cluster"] = pd.concat(scores["Silhouette Score Per Cluster"],axis=1)
        return scores
    
    def results(self,N):
        signatures = self.get_cluster(N)[0]
        if self.pseudosample_on is not None:
            self.sig_vis = SignatureVisualiser.SignatureVisualiser(self.sig_process.data, self.sig_process.metadata,N,True,traditional=self.traditional)
        else:
            self.sig_vis = SignatureVisualiser.SignatureVisualiser(self.sig_process.data, self.sig_process.metadata,N,False,traditional=self.traditional)
            
        scores = self.metrics()
        Results_Table = {}
        for signature in self.sig_vis.signatures:
            Results_Table[str(signature)] = []
            for sig in signatures:
                Results_Table[str(signature)].append(cosine_similarity(X=np.array(self.sig_vis.signatures[signature]).reshape(1, -1),
                                                                       Y=np.array(signatures[sig]).reshape(1, -1))[0][0])
            Results_Table[str(signature)] = [np.argmax(Results_Table[str(signature)]),
                                        Results_Table[str(signature)][np.argmax(Results_Table[str(signature)])],
                                        scores['Silhouette Score Per Cluster'][str(N)][np.argmax(Results_Table[str(signature)])]]
        Results_Table = pd.DataFrame(Results_Table).T
        Results_Table.columns = ["Cluster","Cosine","Silhouette"]
        Results_Table.index.name = "Signature"
        return Results_Table
    
    def get_cluster(self,N):
        catalogue = pd.concat(self.mutational_catalogues[f'{N}'],axis=0)
        exposures = pd.concat(self.mutational_exposures[f'{N}'],axis=0)
        distances = pairwise_distances(catalogue,catalogue,metric=self.cluster_distance)
        clusters = KMeans(n_clusters=N).fit(distances)
        sigs = {}
        expos = {}
        for i,cluster in enumerate(clusters.labels_):
            if str(cluster) not in sigs:
                sigs[str(cluster)] = [catalogue.iloc[i]]
                expos[str(cluster)] = [exposures.iloc[i]]
            else:
                sigs[str(cluster)].append(catalogue.iloc[i])
                expos[str(cluster)].append(exposures.iloc[i])
        signatures = []
        exposures = []
        for key in sigs.keys():
            signatures.append(pd.concat(sigs[key],axis=1).T.mean())
            exposures.append(pd.concat(expos[key],axis=1).T.mean())
        expos = pd.concat(exposures,axis=1)
        signatures = pd.concat(signatures,axis=1)
        expos_absolute = expos
        expos = expos.div(expos.sum(axis=1), axis=0)
        expos.columns = [f'Signature {col}' for col in expos.columns]
        
        ordered_index = pd.Series(expos.index).str.split("_",expand=True)[0]
        expos["order"] = ordered_index.values
        if expos["order"].str.contains("-").any() == True:
            expos["order"] =  pd.to_datetime(expos["order"]).dt.date
        else:
            expos["order"] =  expos["order"].astype('int64')
        expos = expos.groupby(expos.order).mean()
        
        expos_absolute.columns = [f'Signature {col}' for col in expos_absolute.columns]
        expos_absolute = expos_absolute.sort_index()
        return signatures,expos
    
    def get_cluster_signature(self,N):
        return self.get_cluster(N)[0]
    
    def get_cluster_exposure(self,N):
        return self.get_cluster(N)[1]
    
    def plot_exposures(self,exposures,palette=[],figsize=(15,7)):
        fig,ax = plt.subplots(figsize=figsize)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
        if len(palette) == 0:
            exposures.plot(kind="bar",stacked=True,figsize=figsize,width=1,edgecolor="black",title="Signature Exposures per Week",xlabel="Epidemic Week",ylabel="Exposure",ax=ax)
        else:
            exposures.plot(kind="bar",stacked=True,figsize=figsize,width=1,edgecolor="black",title="Signature Exposures per Week",xlabel="Epidemic Week",ylabel="Exposure",ax=ax,color = palette)
        
                
    def plot_signature(self,sig_type,signatures,traditional=False,plot_size=(50,25),comparableY=True):
        cf.plot_signature(sig_type,signatures,traditional,plot_size,comparableY)
        
        
    def plot_mutation_rate(self):
        fig,ax = plt.subplots()
        mutation_rate = pd.concat([self.data,self.metadata["epi_week"]],axis=1)
        mutation_rate["epi_week"] = pd.to_numeric(mutation_rate["epi_week"])
        mutation_rate = mutation_rate.sort_values("epi_week")
        mutation_rate = mutation_rate.groupby("epi_week").sum()
        mutation_rate = mutation_rate.T.groupby([s.split('-')[0] for s in mutation_rate.T.index.values]).sum().T
        palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
        mutation_rate.plot(ax=ax,kind="bar",color=palette,figsize=(30,10),ylabel="Mutations Per Week",stacked=True,width=1.0,title="Mutations per Epidemic Week",xlabel="Epidemic Week",edgecolor='black')
                