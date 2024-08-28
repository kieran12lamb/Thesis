import numpy as np
import pandas as pd
import MutationalSignatures.compiled_functions as cf
from sklearn.decomposition import NMF
import os
import matplotlib.pyplot as plt
import seaborn as sns
import math
import matplotlib.ticker as mtick
from IPython.display import Image
from scipy import stats
from sklearn.metrics.pairwise import cosine_similarity
from pathlib import Path
from Bio.Seq import Seq


class SignatureVisualiser:

    def __init__(self, data, metadata, components,pseudo=False,trimer_proportions=None,traditional=False):
        self.data = data
        self.metadata = metadata
        
        if traditional == True:
            self.make_traditional()
        
        self.decomposition(components)
        
        if trimer_proportions is not None:
            self.normalise_trimers(trimer_proportions)
            
        
        self.signature_probabilities()
        
        self.pseudo = pseudo

    def decomposition(self,components=None, random_state=None):
            import mkl
            mkl.set_num_threads(6)
            #Perform NMF on Matrix
            V = self.data
            #Remove any columns with all zeros from the NMF
            reduced_V = V.loc[:, (V != 0).any(axis=0)]
            nmf = NMF(n_components = components, init="nndsvdar",random_state=random_state,solver="mu",max_iter=5000)
            nmf.fit(reduced_V.T)
            #Perform NMF on Data
            W = nmf.transform(reduced_V.T)
            H = nmf.components_
            self.signatures = pd.DataFrame(W)
            self.signatures.index = reduced_V.columns
            
            zero_V = ((V < 0)*1)[:components].T
            zero_V.columns = self.signatures.columns
            self.signatures = self.signatures.add(zero_V, fill_value=0)
            
            self.exposures = pd.DataFrame(H)
            self.nmf = nmf
        
    def make_traditional(self):
        trad_subs = [ ["CA","CG","CT","TA","TC","TG"],
                        ["GT","GC","GA","AT","AG","AC"] ]
        traditional_sig = []
        data = self.data.T
        for i, sub in enumerate(trad_subs[0]):
            remains = data.iloc[data.index.str.startswith(sub)]
            translate = data.iloc[data.index.str.startswith(trad_subs[1][i])]
            new_index = []
            for index,row in translate.iterrows():
                new_sub =sub
                complementary_codon = str(Seq(index[3:]).complement())
                new_index.append(f"{sub}-{complementary_codon}")
            translate.index = new_index
            translate = translate.sort_index()
            remains = remains.sort_index()
            remains = remains.add(translate)
            traditional_sig.append(remains)
        self.data = pd.concat(traditional_sig,axis=0).T           

    def signature_probabilities(self):
        self.signature_absolute = self.signatures
        self.exposure_absolute = self.exposures
        self.signatures = self.signatures.divide(self.signatures.sum(axis=0)).fillna(0)
        self.exposures = np.array(pd.DataFrame(self.exposures).divide(self.exposures.sum(axis=0)).fillna(0))

    def plot_signature(self,sig_type,signatures,traditional=False,plot_size=(50,25),comparableY=True):
        if sig_type == "Nucleotide":
            max_bar = np.max(np.ravel(signatures.values))
            if traditional:
                trad_subs = [ ["CA","CG","CT","TA","TC","TG"],
                        ["GT","GC","GA","AT","AG","AC"] ]
                traditional_sig = []
                for i, sub in enumerate(trad_subs[0]):
                    traditional_sig.append(np.add(signatures.iloc[signatures.index.str.startswith(sub)],signatures.iloc[signatures.index.str.startswith(trad_subs[1][i])]))
                signatures = pd.concat(traditional_sig,axis=0)
                max_bar = np.max(np.ravel(signatures.values))
                sub_groups = signatures.index.str[:2]
                palette = ['#04BBEC','black','#E42824','grey','#A0CF63','#EEC4C4']
                fig, axes = plt.subplots(len(signatures.columns),len(np.unique(sub_groups)),figsize=(plot_size[0], plot_size[1]), sharey='row')
            else:
                palette = ['#04BBEC','black','#E42824','grey','#A0CF63','#EEC4C4']
                if len(np.unique(signatures.index.str[:2]))>6:
                    palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
                sub_groups = signatures.index.str[:2]
                fig, axes = plt.subplots(len(signatures.columns),len(np.unique(sub_groups)),figsize=(plot_size[0], plot_size[1]), sharey='row')
            print(palette)
            if len(signatures.columns) == 1:
                plot_data = pd.DataFrame(np.array(signatures[signatures.columns[0]]).T)
                plot_data.index = signatures.index
                plot_data.columns = ['values']
                plot_data['sub_groups'] = sub_groups
                plot_data['trimers'] = [trimer[3:] for trimer in plot_data.index]
                for j in range(len(np.unique(sub_groups))):
                    ax = axes[j]
                    ax.set_title(np.unique(sub_groups)[j])
                    sub_plot_data = plot_data[plot_data.sub_groups == np.unique(sub_groups)[j]]
                    sub_plot_data = sub_plot_data.sort_values(by="trimers")
                    sns.barplot(data=sub_plot_data, x='trimers', y="values", hue="sub_groups", palette=[palette[j]],ax=ax,edgecolor='black')
                    ax.get_legend().remove()
                    if comparableY ==False:
                        max_bar = np.max(np.ravel(np.ravel(plot_data["values"])))
                    scaling_factor = max_bar/0.02
                    ax.add_patch(plt.Rectangle((-0.75,max_bar+(0.001*scaling_factor)),25, 0.002*scaling_factor,facecolor=palette[j],clip_on=True,linewidth = 0.1))
                    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
                    ax.set_ylim(0, max_bar+(0.003*scaling_factor))
                    ax.set_xlim(left=-0.75)
                    plt.subplots_adjust(wspace=0,hspace=0.3)
                    ax.set_xticklabels(sub_plot_data.trimers, rotation=90)
                    ax.set_xlabel("")
                    ax.set_ylabel("Mutation Type Probability")
                    ax.set_title(np.unique(sub_groups)[j],fontsize=10)
                    if j !=0:
                        ax.axes.yaxis.set_visible(False)
            else:
                for i in range(len(signatures.columns)):
                    plot_data = pd.DataFrame(np.array(signatures[signatures.columns[i]]).T)
                    plot_data.index = signatures.index
                    plot_data.columns = ['values']
                    plot_data['sub_groups'] = sub_groups
                    plot_data['trimers'] = [trimer[3:] for trimer in plot_data.index]
                    # max_bar = np.max(plot_data['values'])
                    for j in range(len(np.unique(sub_groups))):
                        ax = axes[i][j]
                        ax.set_title(np.unique(sub_groups)[j])
                        sub_plot_data = plot_data[plot_data.sub_groups == np.unique(sub_groups)[j]]
                        sub_plot_data = sub_plot_data.sort_values(by="trimers")
                        sns.barplot(data=sub_plot_data, x='trimers', y="values", hue="sub_groups", palette=[palette[j]],ax=ax,edgecolor='black')
                        ax.get_legend().remove()
                        if comparableY ==False:
                            max_bar = np.max(np.ravel(np.ravel(plot_data["values"])))
                        scaling_factor = max_bar/0.02
                        ax.add_patch(plt.Rectangle((-0.75,max_bar+(0.001*scaling_factor)),25, 0.002*scaling_factor,facecolor=palette[j],clip_on=True,linewidth = 0.1))
                        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
                        ax.set_ylim(0, max_bar+(0.003*scaling_factor))
                        ax.set_xlim(left=-0.75)
                        plt.subplots_adjust(wspace=0,hspace=0.3)
                        ax.set_xticklabels(sub_plot_data.trimers, rotation=90)
                        ax.set_xlabel("")
                        ax.set_ylabel("Mutation Type Probability")
                        ax.set_title(np.unique(sub_groups)[j],fontsize=10)
                        if j !=0:
                            ax.axes.yaxis.set_visible(False)  
                plt.subplots_adjust(wspace=0,hspace=0.3)
        elif sig_type == "Amino Acid":
            cmap = ["#0433FF","#FF2603","#08F900","#000332","#FE35B6","#035301","#FFD300","#009FFF","#9B4D42","#06FBBE","#783FC1","#209698","#FFACFC","#B2CB71","#F1275C","#FE8F42","#DD3DFF","#211A00","#711354","#766C95","#00AD24"]
            max_bar = np.max(np.ravel(signatures.values))
            widths = []
            for i,acid in enumerate(np.sort(signatures.index.str[0].unique())):
                widths.append(len(signatures[signatures.index.str[0] == acid].index.str[2].unique()))
            fig, axes = plt.subplots(len(signatures.columns),len(np.unique(signatures.index.str[0])),figsize=(plot_size[0], plot_size[1]), sharey='row',gridspec_kw={'width_ratios':widths})

            for sig_index in range(len(signatures.columns)):
                plot_data = signatures[signatures.columns[sig_index]]
                summed_vals = plot_data.sum()
                for i,acid in enumerate(np.sort(plot_data.index.str[0].unique())):
                    sub_allAA = plot_data[plot_data.index.str[0] == acid]
                    sub_allAA = sub_allAA.sort_index()
                    axes[sig_index,i].bar(x=sub_allAA.index.str[2], height= sub_allAA.values,color=cmap[i])
                for i, ax in enumerate(axes[sig_index]):
                    ax.set_title(np.sort(plot_data.index.str[0].unique())[i])
                    scaling_factor = max_bar/0.02
                    ax.add_patch(plt.Rectangle((-0.75,max_bar+(0.001*scaling_factor)),
                                                widths[i]*1.005,
                                                0.002*scaling_factor,
                                                facecolor=cmap[i],
                                                linewidth = 0))
                    ax.yaxis.set_major_formatter(mtick.PercentFormatter(summed_vals))
                    ax.set_ylim(0, max_bar+(0.003*scaling_factor))
                    ax.set_xlim(left=-0.75)
                    if i!=0:
                        ax.axes.yaxis.set_visible(False)
                        ax.yaxis.set_major_formatter(mtick.PercentFormatter(summed_vals))
            plt.subplots_adjust(wspace=0,hspace=0.3)
        else:
            print("Sig Type must be Nucleotide or Amino Acid")
    
    def weekly_exposures(self,plot_size=(15,7),comparableY=True):
        fig, axes = plt.subplots(figsize=(plot_size[0], plot_size[1]))
        bar_plot_data = pd.DataFrame()
        no_int = False
        if self.pseudo == False:
            for i,epi_week in enumerate(self.metadata.epi_week.unique()):
                data = pd.DataFrame(self.exposures.T[self.metadata.epi_week == epi_week])
                data = data.mean(axis=0)
                if comparableY:    
                    bar_plot_data[epi_week] = (data/data.sum())*100
                else:
                    bar_plot_data[epi_week] = data
        else:
            if "epi_week" in self.metadata.columns[self.metadata.columns.str.startswith("epi_week")]:
                feature = "epi_week"
            else:
                feature = "sample_date"
                no_int = True
            feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(feature)]
            feature_metadata.columns = [col[len(feature)+1:] for col in feature_metadata.columns]
            if feature in self.metadata.columns:
                for i,time_point in enumerate(self.metadata[feature].unique()):
                    data = pd.DataFrame(self.exposures.T[self.metadata[feature] == time_point])
                    data = data.mean(axis=0)
                    if comparableY:    
                        bar_plot_data[time_point] = (data/data.sum())*100
                    else:
                        bar_plot_data[time_point] = data.copy()
            else:
                for i,epi_week in enumerate(feature_metadata.columns):
                    data = pd.DataFrame(self.exposures.T[feature_metadata[epi_week] >=1])
                    data = data.mean(axis=0)
                    if comparableY:
                        bar_plot_data[epi_week] = (data/data.sum())*100
                    else:
                        bar_plot_data[epi_week] = data.copy()
        bar_plot_data = bar_plot_data.T
        if no_int == False:
            bar_plot_data.index = bar_plot_data.index.astype(int)
        bar_plot_data = bar_plot_data.sort_index().T
        bar_plot_data.index = ["Signature " +str(i) for i in bar_plot_data.index]
#         palette = ["#E05D5D","#00A19D", '#7027A0',"#FFB344",'#F56FAD','#1EAE98',"#111D5E","#FB7813","#E61C5D","#55968F"]
        bar_plot_data.T.plot(kind='bar', stacked=True, ax=axes,width=1,edgecolor='black')
        axes.set_xlabel("Epidemic Week")
        axes.set_ylabel("Percentage of Signature Exposure")
        axes.set_title("Mutation Exposures for each week")
        axes.legend(loc='upper right',bbox_to_anchor=(1.15, 0.5))
        plt.subplots_adjust(wspace=0.3,hspace=0.3)
        axes.yaxis.set_major_formatter(mtick.PercentFormatter())
        cmap = ["#0433FF","#FF2603","#08F900","#000332","#FE35B6","#035301","#FFD300","#009FFF","#9B4D42","#06FBBE","#783FC1","#209698"]
        

        plt.setp(axes.patches, linewidth=1)
        plt.legend()
        plt.show()    

    def export_signatures(self,classes, file_name, traditional=False, exclude_substitution=None,inclusive_substitution=None, signatures = None):
        if signatures is None:
            signatures = self.signatures
            
        if traditional:
            subs = [["CA","CG","CT","TA","TC","TG"],
                    ["GT","GC","GA","AT","AG","AC"]]
            sigs = []
            for i, sub in enumerate(subs[0]):
                sigs.append(np.add(signatures.iloc[signatures.index.str.startswith(sub)],signatures.iloc[signatures.index.str.startswith(subs[1][i])]))
            sigs = pd.concat(sigs,axis=0)
            sigs.columns=["SBS"+str(col)+"_Cov" for col in sigs.columns]
        else:
            sigs = pd.DataFrame(index =self.data.columns, data = signatures).sort_index()
            sigs.columns=["SBS"+str(col)+"_Cov" for col in sigs.columns]
        for signature_name in sigs.columns:
            signature = pd.merge(classes[['sub-trimer','Substitution','Trimer']],sigs[signature_name],right_index=True,left_on="sub-trimer")
            signature = signature.drop("sub-trimer",axis=1)
            signature = signature.rename(columns={"Substitution": "Type", "Trimer": "Subtype"})
            signature['Type'] = [sub[0]+">"+sub[1] for sub in signature['Type']]
            if inclusive_substitution:
                for i in range(len(signature)):
                    print(signature.iloc[i].Type,inclusive_substitution)
                    if signature.iloc[i].Type != inclusive_substitution:
                        signature[signature_name].iloc[i] = 0
            elif exclude_substitution:
                for i in range(len(signature)):
                    if signature.Type.iloc[i] in exclude_substitution:
                        signature[signature_name].iloc[i]=0
            if traditional:
                signature.to_csv(path_or_buf=file_name+signature_name+".csv",index=False)    
            else:
                signature.to_csv(path_or_buf=file_name+signature_name+".csv",index=False)
        
    def exposure_barplot(self,plot_size=(7,7)):
        fig, ax = plt.subplots(figsize=(plot_size[0], plot_size[1]),sharex=True,sharey=True)
        #Normalised weights
        sns.barplot(x=np.arange(len(self.signatures.columns)),y=self.exposures.sum(axis=1)/self.exposures.sum(axis=1).sum(), ax=ax)
        ax.set_ylabel("Attribution Percentage ")
        ax.set_xlabel("Signature")
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1)) 

    def exposure_piecharts(self, interesting_values,feature,plot_size=(15,15)):
        num_plots = len(interesting_values)
        fig, axes = plt.subplots(2,math.ceil(num_plots/2),figsize=(plot_size[0], plot_size[1]),sharex=True,sharey=True)
        axes= axes.flatten()
        if self.pseudo == False:
            stacked = pd.DataFrame(index=self.metadata[feature][self.metadata[feature].isin(interesting_values)],
                                    data=self.exposures.T[self.metadata[feature].isin(interesting_values)])
            stacked = stacked.groupby(level=0).sum()
        else:
            feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(feature)]
            feature_metadata.columns = [col[len(feature)+1:] for col in feature_metadata.columns]
            # print(feature_metadata)
            stacked = pd.DataFrame()
            for i,phenotype in enumerate(feature_metadata.columns):
                    if phenotype in interesting_values:
                        data = pd.DataFrame(self.exposures.T[feature_metadata[phenotype] >=1])
                        data = np.multiply(data,pd.concat([feature_metadata[phenotype][feature_metadata[phenotype] >=1]]*data.shape[1],axis=1))
                        # print(data)
                        data = data.mean(axis=0)
                        stacked[str(phenotype)] = (data/data.sum())*100
            stacked = stacked.T
        for i in range(len(interesting_values)):
            labels = ["Signature "+str(sig) for sig in stacked.columns]
            axes[i].pie(x=stacked.iloc[i], autopct="%.1f%%", explode=[0.15]*len(labels), labels=labels, pctdistance=0.7)
            axes[i].set_title(feature+str(stacked.index[i]), fontsize=14)

    def phenotype_barplots(self, mutation_list,plot_size=(15,15)):

        key_mutations = mutation_list
        if type(mutation_list) is str:
            key_mutation = key_mutations
            fig, axes = plt.subplots(figsize=(plot_size[0], plot_size[1])) 
            bar_plot_data = pd.DataFrame()
            if self.pseudo == False:
                for i,phenotype in enumerate(self.metadata[key_mutation].unique()):
                    data = pd.DataFrame(self.exposures.T[self.metadata[key_mutation] == phenotype])
                    data_mean = data.mean(axis=0)
                    bar_plot_data[str(phenotype)] = (data_mean/data_mean.sum())*100
            else:
                
                feature_metadata = self.metadata.loc[:, np.array(self.metadata.columns.str.startswith(key_mutation))]
                feature_metadata.columns = [col[len(key_mutation)+1:] for col in feature_metadata.columns]
                # print(feature_metadata)
                for i,phenotype in enumerate(feature_metadata.columns):
                    data = pd.DataFrame(self.exposures.T[feature_metadata[phenotype] >=1])
                    data = np.multiply(data,pd.concat([feature_metadata[phenotype][feature_metadata[phenotype] >=1]]*data.shape[1],axis=1))
                    data_mean = data.mean(axis=0)
                    bar_plot_data[str(phenotype)] = (data_mean/data_mean.sum())*100
            bar_plot_data.index = ["Signature " +str(i) for i in bar_plot_data.index]
            bar_plot_data.T.plot(kind='bar', stacked=True, ax=axes)
            axes.set_xlabel("Amino Acid Phenotype")
            axes.set_ylabel("Percentage of Signature Exposure")
            axes.set_title("Mutation Exposures for each Phenotype of "+key_mutation)
            axes.legend(loc='upper right',bbox_to_anchor=(1.15, 0.5))
            plt.subplots_adjust(wspace=0.3,hspace=0.3)
            axes.yaxis.set_major_formatter(mtick.PercentFormatter()) 
        else:
            fig, axes = plt.subplots(math.ceil(len(key_mutations)/2),2,figsize=(plot_size[0], plot_size[1]))
            axes = axes.flatten()
            for ax_index, key_mutation in enumerate(key_mutations):
                bar_plot_data = pd.DataFrame()
                if self.pseudo == False:
                    for i,phenotype in enumerate(self.metadata[key_mutation].unique()):
                        data = pd.DataFrame(self.exposures.T[self.metadata[key_mutation] == phenotype])
                        data_mean = data.mean(axis=0)
                        bar_plot_data[str(phenotype)] = (data_mean/data_mean.sum())*100
                else:
                    feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(key_mutation)]
                    feature_metadata.columns = [col[len(key_mutation)+1:] for col in feature_metadata.columns]
                    # print(feature_metadata)
                    for i,phenotype in enumerate(feature_metadata.columns):
                        data = pd.DataFrame(self.exposures.T[feature_metadata[phenotype] >=1])
                        data = np.multiply(data,pd.concat([feature_metadata[phenotype][feature_metadata[phenotype] >=1]]*data.shape[1],axis=1))
                        data_mean = data.mean(axis=0)
                        bar_plot_data[str(phenotype)] = (data_mean/data_mean.sum())*100
                bar_plot_data.index = ["Signature " +str(i) for i in bar_plot_data.index]
                bar_plot_data.T.plot(kind='bar', stacked=True, ax=axes[ax_index])
                axes[ax_index].set_xlabel("Amino Acid Phenotype")
                axes[ax_index].set_ylabel("Percentage of Signature Exposure")
                axes[ax_index].set_title("Mutation Exposures for each Phenotype of "+key_mutation)
                axes[ax_index].legend(loc='upper right',bbox_to_anchor=(1.15, 0.5))
                plt.subplots_adjust(wspace=0.3,hspace=0.3)
                axes[ax_index].yaxis.set_major_formatter(mtick.PercentFormatter()) 

    def exposure_boxplots(self,mutation,plot_size=(30,15)):
        if self.pseudo == False:
            num_plots = len(np.unique(self.metadata[mutation]))
            mutation = self.metadata[mutation]
            phenotypes = np.unique(mutation)
        else:
            feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(mutation)]
            feature_metadata.columns = [col[len(mutation)+1:] for col in feature_metadata.columns]
            num_plots = len(feature_metadata.columns)
            phenotypes = feature_metadata.columns
        fig, axes = plt.subplots(2,math.ceil(num_plots/2),figsize=(plot_size[0], plot_size[1]))
        axes = axes.flatten()

        for ax_index, phenotype in enumerate(phenotypes):
            if self.pseudo == False:
                data = pd.DataFrame(self.exposures.T[mutation == phenotype])
            else:
                data = pd.DataFrame(self.exposures.T[feature_metadata[phenotype] >=1 ])
            sns.boxplot(data=data,ax=axes[ax_index])
            axes[ax_index].set_xlabel("Signature")
            axes[ax_index].set_ylabel("Exposure")
            axes[ax_index].set_title(phenotype)
            axes[ax_index].yaxis.set_major_formatter(mtick.PercentFormatter(1))

    def sig2aa_plot(self,signature_file,output,genome_type):
        os.system('Rscript --vanilla sigtoaa.R ' + signature_file +" "+output+" "+genome_type)
        return Image(filename=(output))

    def signature_orf_risk(self,reference_sequence,orfs,signature,normalise_orfs=False,plot_size=(15,5)):
        signature = pd.DataFrame(signature[signature!=0])
        signature.columns = ["proportion"]
        signature["substitution"] = signature.index.str.slice(stop=2)
        signature["trimer"] = signature.index.str.slice(start=3)
        signature_trimers = np.unique(signature.index.str.slice(start=3))
        fig, axes = plt.subplots(figsize=(plot_size[0], plot_size[1]))
        orf_risks = []

        for i, orf in enumerate(orfs.index):
            orf_risk = []
            orf_nucleotides = orfs.iloc[i].Locations.extract(reference_sequence)
            for substitution in signature.substitution.unique():
                substitution_filtered_signature = signature[signature.substitution == substitution]
                for trimer in substitution_filtered_signature.trimer.unique():
                    signature_proportion = substitution_filtered_signature[substitution_filtered_signature.trimer ==trimer]["proportion"].sum()
                    if normalise_orfs == True:
                        trimer_proportion = np.divide(np.count_nonzero([orf_nucleotides[base_index:base_index+3] == trimer for base_index in range(len(orf_nucleotides)-2)]),(len(orf_nucleotides)-2))
                        total_risk = trimer_proportion*signature_proportion
                    else:
                        trimer_count = np.count_nonzero([orf_nucleotides[base_index:base_index+3] == trimer for base_index in range(len(orf_nucleotides)-2)])
                        total_risk = trimer_count*signature_proportion
                    orf_risk.append(pd.DataFrame({f'{substitution}':total_risk},index=[0]))
                    
            orf_risk = pd.concat(orf_risk,axis=1).T
            orf_risk.columns = ["ORF-"+str(i)+": "+orf] 
            orf_risks.append(orf_risk)
        orf_risks = pd.concat(orf_risks,axis=1)
        orf_risks = orf_risks.T.sum(axis=1, level=0)
        orf_risks = orf_risks.T.sort_index().T
        palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
        orf_risks.plot(kind="bar",ax=axes,stacked=True,color = palette)
        # axes.yaxis.set_major_formatter(mtick.PercentFormatter(1))
        axes.set_xlabel("Trimers")
        axes.set_ylabel("ORF Enrichment")
        axes.legend(loc='upper right',bbox_to_anchor=(1.15, 0.5))
        plt.xticks(rotation=0)
        
    def plot_mutation_rate_per_signature(self,figsize=(30,10)):
        exposures = (pd.DataFrame(self.exposures)/pd.DataFrame(self.exposures).sum())*100
        mutation_rate = pd.concat([self.data,self.metadata["epi_week"]],axis=1)
        mutation_rate["epi_week"] = pd.to_numeric(mutation_rate["epi_week"])
        mutation_rate = mutation_rate.sort_values("epi_week")
        mutation_rate = mutation_rate.groupby("epi_week").sum().sum(axis=1)
#         mutation_rate = mutation_rate.T.groupby([s.split('-')[0] for s in mutation_rate.T.index.values]).sum().T

        palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
        
        fig, axes = plt.subplots(figsize=figsize)
        bar_plot_data = pd.DataFrame()
        no_int = False
        if self.pseudo == False:
            for i,epi_week in enumerate(self.metadata.epi_week.unique()):
                data = pd.DataFrame(self.exposures.T[self.metadata.epi_week == epi_week])
                data = data.mean(axis=0)
                bar_plot_data[epi_week] = data
        else:
            if len(self.metadata.columns.str.startswith("epi_week") == True)>1:
                feature = "epi_week"
            else:
                feature = "sample_date"
                no_int = True
            feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(feature)]
            feature_metadata.columns = [col[len(feature)+1:] for col in feature_metadata.columns]
            if feature in self.metadata.columns:
                for i,time_point in enumerate(self.metadata[feature].unique()):
                    data = pd.DataFrame(self.exposures.T[self.metadata[feature] == time_point])
                    data = data.mean(axis=0)
                    bar_plot_data[time_point] = data.copy()
            else:
                for i,epi_week in enumerate(feature_metadata.columns):
                    data = pd.DataFrame(self.exposures.T[feature_metadata[epi_week] >=1])
                    data = data.mean(axis=0)
                    bar_plot_data[epi_week] = data.copy()
        bar_plot_data = bar_plot_data.T
        if no_int == False:
            bar_plot_data.index = bar_plot_data.index.astype(int)
        bar_plot_data = bar_plot_data.sort_index().T
        bar_plot_data.index = [str(i) for i in bar_plot_data.index]
#         palette = ["#E05D5D","#00A19D", '#7027A0',"#FFB344",'#F56FAD','#1EAE98',"#111D5E","#FB7813","#E61C5D","#55968F"]
        
        cmap = ["#0433FF","#FF2603","#08F900","#000332","#FE35B6","#035301","#FFD300","#009FFF","#9B4D42","#06FBBE","#783FC1","#209698"]
        sig_exposure_count_plot = (mutation_rate*bar_plot_data).T
        sig_exposure_count_plot.columns = [f'Signature {col}'for col in sig_exposure_count_plot.columns]
        sig_exposure_count_plot.plot(kind="bar", stacked=True,ax=axes,width=1.0,edgecolor='black') 
        axes.set_xlabel("Epidemic Week")
        axes.set_ylabel("Mutations")
        axes.set_title("Mutation Count per Signature")
        plt.subplots_adjust(wspace=0.3,hspace=0.3)
        
        
    def plot_mutation_rate(self):
        fig,ax = plt.subplots()
        mutation_rate = pd.concat([self.data,self.metadata["epi_week"]],axis=1)
        mutation_rate["epi_week"] = pd.to_numeric(mutation_rate["epi_week"])
        mutation_rate = mutation_rate.sort_values("epi_week")
        mutation_rate = mutation_rate.groupby("epi_week").sum()
        mutation_rate = mutation_rate.T.groupby([s.split('-')[0] for s in mutation_rate.T.index.values]).sum().T
        palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
        print(mutation_rate)
        mutation_rate.plot(ax=ax,kind="bar",color=palette,figsize=(30,10),ylabel="Mutations Per Week",stacked=True,width=1.0,title="Mutations per Epidemic Week",xlabel="Epidemic Week",edgecolor='black')
        
    def export_to_observable(self,comparableY=True):
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, r'Resources/Observable_Exports')
        if not os.path.exists(final_directory):
           os.makedirs(final_directory)
        #Export Signatures
        sig_export = ((self.signatures/self.signatures.sum())*100)
        sig_export.index.name = "Mutation"
        labels = pd.DataFrame(pd.Series(sig_export.index).str.split("-",expand = True))
        labels.index = sig_export.index
        labels.columns = ["SUB","CONTEXT"]
        sig_export.columns = [f"Signature_{sig_num}"for sig_num in sig_export]
        sig_export = sig_export/100
        sig_export = pd.merge(labels,sig_export,how="left",left_index=True,right_index=True)
        sig_export = sig_export.sort_index()
        print(sig_export)
        sig_export.to_csv("Resources/Observable_Exports/Signatures.csv")
        # Export Exposures
        bar_plot_data = pd.DataFrame()
        no_int = False
        if self.pseudo == False:
            for i,epi_week in enumerate(self.metadata.epi_week.unique()):
                data = pd.DataFrame(self.exposures.T[self.metadata.epi_week == epi_week])
                data = data.mean(axis=0)
                if comparableY:    
                    bar_plot_data[epi_week] = (data/data.sum())*100
                else:
                    bar_plot_data[epi_week] = data
        else:
            if len(self.metadata.columns.str.startswith("epi_week") == True)>1:
                feature = "epi_week"
            else:
                feature = "sample_date"
                no_int = True
            feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(feature)]
            feature_metadata.columns = [col[len(feature)+1:] for col in feature_metadata.columns]
            if feature in self.metadata.columns:
                for i,time_point in enumerate(self.metadata[feature].unique()):
                    data = pd.DataFrame(self.exposures.T[self.metadata[feature] == time_point])
                    data = data.mean(axis=0)
                    if comparableY:    
                        bar_plot_data[time_point] = (data/data.sum())*100
                    else:
                        bar_plot_data[time_point] = data.copy()
            else:
                for i,epi_week in enumerate(feature_metadata.columns):
                    data = pd.DataFrame(self.exposures.T[feature_metadata[epi_week] >=1])
                    data = data.mean(axis=0)
                    if comparableY:
                        bar_plot_data[epi_week] = (data/data.sum())*100
                    else:
                        bar_plot_data[epi_week] = data.copy()
        bar_plot_data = bar_plot_data.T
        if no_int == False:
            bar_plot_data.index = bar_plot_data.index.astype(int)
        bar_plot_data = bar_plot_data.sort_index().T
        bar_plot_data.index = ["Signature " +str(i) for i in bar_plot_data.index]
#         palette = ["#E05D5D","#00A19D", '#7027A0',"#FFB344",'#F56FAD','#1EAE98',"#111D5E","#FB7813","#E61C5D","#55968F"]
        
        exposure_export = bar_plot_data.T.melt(value_vars=bar_plot_data.T.columns,ignore_index=False)
        if feature == "epi_week":
            feature = "Epi_Week"
        exposure_export.index.name = feature
        exposure_export.columns = ["Signature","Exposure"]
        exposure_export["Exposure"]= exposure_export["Exposure"]/100
        print(exposure_export)
        exposure_export.to_csv("Resources/Observable_Exports/Exposures.csv")
        
        
        
        exposures = (pd.DataFrame(self.exposures)/pd.DataFrame(self.exposures).sum())*100
        mutation_rate = pd.concat([self.data,self.metadata["epi_week"]],axis=1)
        mutation_rate["epi_week"] = pd.to_numeric(mutation_rate["epi_week"])
        mutation_rate = mutation_rate.sort_values("epi_week")
        mutation_rate = mutation_rate.groupby("epi_week").sum().sum(axis=1)
        print(mutation_rate)
#         mutation_rate = mutation_rate.T.groupby([s.split('-')[0] for s in mutation_rate.T.index.values]).sum().T

        palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
        print((pd.DataFrame(exposures)/pd.DataFrame(exposures).sum()).T)
        
        bar_plot_data = pd.DataFrame()
        no_int = False
        if self.pseudo == False:
            for i,epi_week in enumerate(self.metadata.epi_week.unique()):
                data = pd.DataFrame(self.exposures.T[self.metadata.epi_week == epi_week])
                data = data.mean(axis=0)
                bar_plot_data[epi_week] = data
        else:
            if len(self.metadata.columns.str.startswith("epi_week") == True)>1:
                feature = "epi_week"
            else:
                feature = "sample_date"
                no_int = True
            feature_metadata = self.metadata.loc[:, self.metadata.columns.str.startswith(feature)]
            feature_metadata.columns = [col[len(feature)+1:] for col in feature_metadata.columns]
            if feature in self.metadata.columns:
                for i,time_point in enumerate(self.metadata[feature].unique()):
                    data = pd.DataFrame(self.exposures.T[self.metadata[feature] == time_point])
                    data = data.mean(axis=0)
                    bar_plot_data[time_point] = data.copy()
            else:
                for i,epi_week in enumerate(feature_metadata.columns):
                    data = pd.DataFrame(self.exposures.T[feature_metadata[epi_week] >=1])
                    data = data.mean(axis=0)
                    bar_plot_data[epi_week] = data.copy()
        bar_plot_data = bar_plot_data.T
        if no_int == False:
            bar_plot_data.index = bar_plot_data.index.astype(int)
        bar_plot_data = bar_plot_data.sort_index().T
        bar_plot_data.index = [str(i) for i in bar_plot_data.index]
        sig_exposure_count_plot = (mutation_rate*bar_plot_data).T
        sig_exposure_count_plot.columns = [f'Signature {col}'for col in sig_exposure_count_plot.columns]
        sig_exposure_count_plot.index.name = "Epi_Week"
        sig_exposure_count_plot = sig_exposure_count_plot.melt(value_vars=sig_exposure_count_plot.columns,ignore_index=False)
        sig_exposure_count_plot.index.name = "Epi_Week"
        sig_exposure_count_plot.columns = ["Signature","Exposure"]
        print(sig_exposure_count_plot)
        sig_exposure_count_plot.to_csv("Resources/Observable_Exports/Exposure_Count.csv")
        
        
        
        mutation_rate = pd.concat([self.data,self.metadata["epi_week"]],axis=1)
        mutation_rate["epi_week"] = pd.to_numeric(mutation_rate["epi_week"])
        mutation_rate = mutation_rate.sort_values("epi_week")
        mutation_rate = mutation_rate.groupby("epi_week").sum()
        mutation_rate = mutation_rate.T.groupby([s.split('-')[0] for s in mutation_rate.T.index.values]).sum().T
        mutation_rate = mutation_rate.melt(value_vars=mutation_rate.columns,ignore_index=False)
        mutation_rate.columns = ["Substitution","Count"]
        mutation_rate.index.name = "Epi_Week"
        mutation_rate.to_csv("Resources/Observable_Exports/Substitution_Counts.csv")
        print(mutation_rate)