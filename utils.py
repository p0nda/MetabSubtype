from scipy.stats import mannwhitneyu,wilcoxon,shapiro, levene, ttest_ind,f_oneway,kruskal
import pingouin as pg
from itertools import combinations
import lifelines
import numpy as np
import pandas as pd
import ast
from collections import *
import pickle
from xgboost import XGBClassifier, XGBRegressor
import xgboost
import shap
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, normalize, StandardScaler
from sklearn.metrics import roc_curve, auc, accuracy_score, make_scorer
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from typing import List
import warnings
import numpy as np
import pprint
from abc import ABC, abstractmethod
warnings.filterwarnings("ignore")
shap.initjs()

# DRAW
def draw_roc(df,selected_features,class_laebl,FOLD_NUM=5,dirpath=''):
    df=df[~(df[class_label].isna())]
    X, y = get_X_and_y(df,[class_label],[])
    features = selected_features
    # feature selection
    X = X.loc[:, features]
    xgb_model = XGB('XGB', {})
    if FOLD_NUM<1:
        plt.figure()
        SEED=13
        X_train, X_test, y_train, y_test = train_test_split(
            X,y,test_size=FOLD_NUM,random_state=SEED)
        xgb_model = XGB('XGB', {})
        xgb_model.train_model(X_train,y_train)
        pred_proba = xgb_model.model.predict_proba(X_test)
        fpr, tpr, _ = roc_curve(y_test, pred_proba[::,1])
        auc_value = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=2,
                    label='ROC curve - (area = %0.4f)' % auc_value)
        # aucs = [element for lis in aucs for element in lis]
        # print(fprs)
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([-0.01, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'response: {features}')
        plt.legend(loc="lower right")
        # print(pred_proba)
    else:
        xgb_model.run_k_fold(X, y,n_splits=FOLD_NUM)
        for fold_index in range(FOLD_NUM):
            fprs = xgb_model.cv_storage.fprs[fold_index]
            tprs = xgb_model.cv_storage.tprs[fold_index]
            aucs = xgb_model.cv_storage.aucs[fold_index]
            auc_value = np.mean(aucs)
            plt.plot(fprs, tprs, lw=2,
                        label=str(fold_index) + '_ROC curve - (area = %0.4f)' % auc_value)
            # aucs = [element for lis in aucs for element in lis]
            # print(fprs)
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(features)
        plt.legend(loc="lower right")
    if dirpath=='':
        plt.show()
    else :
        features=str(features).replace(':','_')
        features=str(features).replace('.','_')
        filepath=dirpath+features
        plt.savefig(filepath)

def draw_box(df,target_features,class_label,dirpath=''):
    # determine boundary
    basic_cols=[class_label]
    df= df.sort_values(class_label)
    # display(df[class_label])
    df.index=[i for i in range(df.shape[0])]
    count_dict=dict(df.value_counts(class_label))
    count_dict=dict(sorted(count_dict.items()))
    class_num=len(count_dict.keys())
    BOUND_INDEXES=[0 for i in range(class_num)]
    for i in range(class_num):
        count_values=list(count_dict.values())
        BOUND_INDEXES[i]=0
        if i==0:
            BOUND_INDEXES[i]=count_values[i]
        else:
            for j in range(i+1):
                BOUND_INDEXES[i]+=count_values[j]
    BOUND_INDEXES.insert(0,0)
    pvalue_dict=get_pvalue_dict(df,[class_label])
    for index,feature in enumerate(target_features):
        # plt.rcParams.update({'font.size': 10})
        plt.figure(figsize=(3,4))
        feature_name=feature
        plt.title(feature_name,fontsize=12)
        target_cols=[]
        target_cols.extend(basic_cols)
        target_cols.append(feature)
        target_df=df.loc[:,target_cols]
        # core step to make sure sorting process is right
        target_df[[feature]]=target_df[[feature]].astype(float)
        # target_df=target_df.groupby(class_label).apply(lambda x: x.sort_values(target_df.columns[1], ascending=True))
        target_df=target_df.sort_values(by=basic_cols)
        # display(target_df)
        y=np.array(target_df[feature].values.astype(float))
        # print(f'{feature}:{y-target_df[feature].values}')
        ## normalize
        y=(y-np.min(y))/(np.max(y)-np.min(y))
        grouped_data=[]
        for i in range(class_num):
            grouped_data.append(y[BOUND_INDEXES[i]:BOUND_INDEXES[i+1]])

        # draw
        labels=list(count_dict.keys())
        if ('T' in labels): 
            if (labels.index('T')!=0):
                tmp_data=grouped_data[-1]
                grouped_data[-1]=grouped_data[0]
                grouped_data[0]=tmp_data
                labels=list(reversed(labels))
        # plt.rcParams.update({'axes.labelsize':19})
        for i in range(class_num):
            plt.scatter(len(grouped_data[i])*[i+1],grouped_data[i])
        plt.boxplot(grouped_data,widths=0.6,labels=labels)
        p_value=pvalue_dict[feature]
        # p_value=0
        plt.legend(["p=%.4f" % p_value],loc='best',handlelength=0,handletextpad=0)
        # save picture
        if dirpath=='':
            plt.show()  
        else:
            filepath=dirpath+feature
            filepath=filepath.replace(':','-')
            plt.savefig(filepath)

def draw_box_with_test(df,target_features,class_label,order=None,color_map=None):
    pvalue_dict=get_pvalue_dict(df,[class_label])
    if (color_map==None) | (order==None):
        using_colors=None
    else:
        using_colors=[color_map[i] for i in order]

    for index,feature_name in enumerate(target_features):
        fig,ax=plt.subplots()
        fig.set_size_inches(4,4)
        ax.set_title(feature_name,fontsize=12)
        # ax=sns.boxplot(df[class_label],df[feature],order=order,color='white',medianprops={"color": "coral"},)
        sns.boxplot(data=df,x=class_label,y=feature_name,order=order,width=0.6,color='white',ax=ax)

        statannot.add_stat_annotation(ax, data=df,x=feature_name,y=class_label ,order=order,
                            box_pairs=[("I", "A"), ("I", "C"), ("I", "N")],
                            test='Mann-Whitney', text_format='star', verbose=2)
        ax.set_ylabel('log2')
        sns.scatterplot(data=df,x=class_label,y=feature_name,hue=class_label,s=50,palette=using_colors,legend=False,ax=ax)
        p_value=pvalue_dict[feature_name]
        ax.legend(["p=%.4f" % p_value], loc='best', borderaxespad=0,handlelength=0,handletextpad=0)

        if dirpath=='':
            plt.show()  
        else:
            feature_name=feature_name.replace(':','-')
            feature_name=feature_name.replace('.','-')
            filepath=dirpath+feature_name
            plt.savefig(filepath)
# get p_value
def get_col_pvalue(df,target_cols, label_cols=['cure'],unwanted_cols=[]):
    # delete unwanted
    pvalue_dict={}
    classes=list(dict(df.value_counts(label_cols[0])).keys())
    colnames=list(set(target_cols)&set(df.columns))
    for i, colname in enumerate(colnames):

        feature_0 = df.loc[(df[label_cols[0]] == classes[0])].loc[:, colname].tolist()
        feature_1 = df.loc[(df[label_cols[0]] == classes[1])].loc[:, colname].tolist()
        feature_1 = list(map(float, feature_1))
        feature_0 = list(map(float, feature_0))
        # display mode
        # print('No.{0} feature: {1}'.format(i,colname) )
        # print(feature_1)
        # print(feature_0)
        is_norm=True
        is_equal_var=True
        if (shapiro(feature_0).pvalue < 0.05 or shapiro(feature_1).pvalue < 0.05):
            cols_not_norm.append(colname)
            is_norm=False
        if (levene(feature_1, feature_0).pvalue > 0.05):
            is_equal_var = True
        else:
            is_equal_var = False
        if is_norm:
            p_value=ttest_ind(feature_1, feature_0, equal_var=is_equal_var).pvalue 
        if is_norm==False:  
            p_value=mannwhitneyu(feature_0,feature_1).pvalue
        pvalue_dict[colname]=p_value
    return pvalue_dict

# 二分类
def get_important_df_ttest(df, label_cols=['cure'],unwanted_cols=[]):
    # delete unwanted
    colnames = df.columns
    colnames = list(set(colnames)-set(label_cols)-set(unwanted_cols))
    classes=list(dict(df.value_counts(label_cols[0])).keys())
    cols_not_important = []
    cols_important = []
    cols_not_norm = []
    for i, colname in enumerate(colnames):

        feature_0 = df.loc[(df[label_cols[0]] == classes[0])].loc[:, colname].tolist()
        feature_1 = df.loc[(df[label_cols[0]] == classes[1])].loc[:, colname].tolist()
        feature_1 = list(map(float, feature_1))
        feature_0 = list(map(float, feature_0))
        # display mode
        # print('No.{0} feature: {1}'.format(i,colname) )
        # print(feature_1)
        # print(feature_0)
        is_norm=True
        is_equal_var=True
        if (shapiro(feature_0).pvalue < 0.05 or shapiro(feature_1).pvalue < 0.05):
            cols_not_norm.append(colname)
            is_norm=False
        if (levene(feature_1, feature_0).pvalue > 0.05):
            is_equal_var = True
        else:
            is_equal_var = False
        if is_norm:
            if ttest_ind(feature_1, feature_0, equal_var=is_equal_var).pvalue < 0.05:
                cols_important.append(colname)
            else:
                cols_not_important.append(colname)
        if is_norm==False:  
            p_value=mannwhitneyu(feature_0,feature_1).pvalue
            if p_value>0.05:
                cols_not_important.append(colname)

    df = df.drop(columns=cols_not_important)
    return df

# 二分类
def get_pvalue_df(df, label_cols=['cure'],unwanted_cols=[]):
    # delete unwanted
    colnames = df.columns
    colnames = list(set(colnames)-set(label_cols)-set(unwanted_cols))
    classes=list(dict(df.value_counts(label_cols[0])).keys())
    cols_not_important = []
    cols_important = []
    cols_not_norm = []
    pvalue_dict={}
    for i, colname in enumerate(colnames):

        feature_0 = df.loc[(df[label_cols[0]] == classes[0])].loc[:, colname].tolist()
        feature_1 = df.loc[(df[label_cols[0]] == classes[1])].loc[:, colname].tolist()
        feature_1 = list(map(float, feature_1))
        feature_0 = list(map(float, feature_0))
        # display mode
        # print('No.{0} feature: {1}'.format(i,colname) )
        # print(feature_1)
        # print(feature_0)
        is_norm=True
        is_equal_var=True
        if (shapiro(feature_0).pvalue < 0.05 or shapiro(feature_1).pvalue < 0.05):
            cols_not_norm.append(colname)
            is_norm=False
        if (levene(feature_1, feature_0).pvalue > 0.05):
            is_equal_var = True
        else:
            is_equal_var = False
        if is_norm:
            col_pvalue=ttest_ind(feature_1, feature_0, equal_var=is_equal_var).pvalue
            if  col_pvalue< 0.05:
                cols_important.append(colname)
            else:
                cols_not_important.append(colname)
        if is_norm==False:  
            col_pvalue=mannwhitneyu(feature_0,feature_1).pvalue
            if col_pvalue>0.05:
                cols_not_important.append(colname)
        pvalue_dict[colname]=col_pvalue
    df=df.append(pvalue_dict,ignore_index=True)
    # df = df.drop(columns=cols_not_important)
    return df
# 多分类
def get_significant_df(df,label_cols=[],unwanted_cols=['id','cure']):
    print(unwanted_cols)
    # delete unwanted
    colnames = df.columns
    colnames = list(set(colnames)-set(label_cols)-set(unwanted_cols))
    classes=list(dict(df.value_counts(label_cols[0])).keys())
    print(classes)
    cols_not_important = []
    cols_important = []
    cols_not_norm = []
    for i,colname in enumerate(colnames):
        class_features=[]
        for label in classes:
            tmp_feature=[]
            # tmp_feature = df.loc[(df[label_cols[0]] == label)]
            tmp_feature = df.loc[(df[label_cols[0]] == label)].loc[:, colname].tolist()
            class_features.append(tmp_feature)
        is_norm=True
        is_equal_var=True
        # data type convert
        df=df.astype({colname:float})
        # judge norm and equal variance
        for class_feature in class_features:
            if shapiro(class_feature).pvalue<0.05:
                is_norm=False
                break
        if levene(*class_features).pvalue<0.05:
            is_equal_var=False
        # select functions by norm and var equality
        if is_equal_var and is_norm:
            _, p_value = f_oneway(*class_features)
            if p_value > 0.05:
                cols_not_important.append(colname)
        elif is_norm and (not is_equal_var):
            p_value=pg.welch_anova(dv=colname,between=label_cols[0],data=df)['p-unc'].values[0]
            if p_value>0.05:
                cols_not_important.append(colname)
        else:
            _, p_value = kruskal(*class_features)
            if p_value>0.05:
                cols_not_important.append(colname)
        cols_important = list(set(df.columns)-set(cols_not_important))
        df = df.loc[:, cols_important]
    return df

# 二分类和多分类

def get_pvalue_dict(df,label_cols,unwanted_cols=[]):
    
    # delete unwanted
    colnames = df.columns
    colnames = list(set(colnames)-set(label_cols)-set(unwanted_cols))
    classes=list(dict(sorted(dict(df.value_counts(label_cols[0])).items())).keys())
    print(classes)
    pvalue_dict={}
    if(len(classes)>2):
        for i,colname in enumerate(colnames):
            class_features=[]
            df=df.astype({colname:float})
            for label in classes:
                tmp_feature=[]
                # tmp_feature = df.loc[(df[label_cols[0]] == label)]
                tmp_feature = df.loc[(df[label_cols[0]] == label)].loc[:, colname].tolist()
                class_features.append(tmp_feature)

            _, p_value = kruskal(*class_features)
            pvalue_dict[colname]=p_value
    else:
        for i,colname in enumerate(colnames):
            class_features=[]
            df=df.astype({colname:float})
            for label in classes:
                tmp_feature=[]
                # tmp_feature = df.loc[(df[label_cols[0]] == label)]
                tmp_feature = df.loc[(df[label_cols[0]] == label)].loc[:, colname].tolist()
                class_features.append(tmp_feature)

            _, p_value = mannwhitneyu(*class_features)
            pvalue_dict[colname]=p_value
        
    return pvalue_dict

# make folder
def mkdir(path):
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径
        print ("---  new folder...  ---")
        print ("---  OK  ---")
    else:
        print ("---  There is this folder!  ---")
# X,y 
def get_X_and_y(df,label_cols=[],unwanted_cols=[]):
    # mapping str to num
    unwanted_cols.extend(label_cols)
    X=df.drop(columns=unwanted_cols,axis=1)
    X=X.astype('float')
    y=df[label_cols[0]].values
    mapping=dict(zip(sorted(set(y)),range(len(y))))
    y=[mapping[i] for i in y]
    y=np.array(y)    
    return X,y
# normalize 
def normalize_X(df,unwanted_cols=[]):
    target_cols=set(df.columns)-set(unwanted_cols)
    X=df.loc[:,target_cols]
    X=X.astype(float)
    X=(X-X.mean())/(X.std())
    df=pd.merge(X,df.loc[:,unwanted_cols],left_index=True,right_index=True)
    return df
# display utils
def display_all():
    pd.set_option('display.max_rows',None)
    pd.set_option('display.max_columns',None)
def display_part():
    pd.set_option('display.max_rows',8)
    pd.set_option('display.max_columns',10)

# analyze combination result

def get_feature_sets(df,threshold=0.95,feature_col='feature'):
    feature_sets=[]
    target_df=df[(df.auc>=threshold)&(df.accuracy>=threshold)]
    tmp_features=target_df[feature_col].values.tolist()
    for features in tmp_features:
        features=ast.literal_eval(features)
        feature_sets.append(features)
    return feature_sets
def get_feature_counter(df,threshold=0.95,feature_col='feature'):
    all_features=[]
    target_df=df[(df.auc>=threshold)&(df.accuracy>=threshold)]
    tmp_features=target_df[feature_col].values.tolist()
    for features in tmp_features:
        features=ast.literal_eval(features)
        all_features.extend(features)
    all_feature_counter=Counter(all_features)
    all_feature_counter=dict(all_feature_counter)
    all_feature_counter=dict(sorted(all_feature_counter.items(),key=lambda item:item[1],reverse=True))
    return all_feature_counter
    all_features=list(all_feature_counter.keys())
    len(set(all_features))

def find_differential_metabolites(df, metabolite_num,class_label,batch_cluster=0):
    """
    Finds differential metabolites of a batch of samples based on K-means clustering.
    
    Args:
    - df: a pandas DataFrame in shape (50,200) storing 50 samples, the last column of the df `group` is used to label the classification of samples.
    - class_label: a string representing the label of the batch of samples
    
    Returns:
    - a list of differential metabolites
    """
    
    # Set the number of clusters to the number of unique labels in the class_label column
    k = len(df[class_label].unique())
    
    # # Perform K-means clustering
    # kmeans = KMeans(n_clusters=k, random_state=0).fit(df.iloc[:, :-1])
    
    # # Find the samples that belong to the batch of interest
    # batch_samples = df[df[class_label] == class_label].index
    
    # # Find the cluster that contains the most batch samples
    # batch_cluster = kmeans.predict(df.loc[batch_samples, :-1].mean().values.reshape(1, -1))[0]
    
    # Find the samples that belong to the batch cluster
    batch_cluster_samples = df[df[class_label] == batch_cluster].index
    
    # Find the differential metabolites between the batch cluster and the rest of the samples
    differential_metabolites = []
    for metabolite in df.columns[:metabolite_num]:
        batch_cluster_metabolite = df.loc[batch_cluster_samples, metabolite]
        other_cluster_metabolite = df.loc[~df.index.isin(batch_cluster_samples), metabolite]
        p_value = ttest_ind(batch_cluster_metabolite, other_cluster_metabolite)[1]
        if p_value < 0.05:
            differential_metabolites.append(metabolite)
    
    return differential_metabolites

# Outliner Process


def get_IQR(x):
    return (np.quantile(x, q=0.75) - np.quantile(x, q=0.25))

def count_outliers(x,iqr_time=1.5):
    outlier_num = len(x[x > (np.quantile(x, q=0.75) + (get_IQR(x)) * iqr_time)])
    return outlier_num

def replace_outliers(x,iqr_time=1.5):
    threshold = (np.quantile(x, q=0.75) + (get_IQR(x)) * iqr_time)
    x[x > threshold] = np.mean(x[x < threshold])
    return x

def delete_outliner(df,metab_num,iqr_time=1.5):
    rawdata_df = df.iloc[:, :metab_num]
    nums = rawdata_df.apply(lambda x : count_outliers(x,iqr_time), axis=0)
    # print(nums)
    fig=plt.figure()
    fig.add_subplot(121)
    plt.hist(nums)
    abnormal_metabolites=[]
    for index,num in enumerate(nums):
        if num>3:
            abnormal_metabolites.append(rawdata_df.columns[index])
    delete_outliner_df=rawdata_df.drop(columns=abnormal_metabolites)
    b = delete_outliner_df.apply(lambda x :replace_outliers(x,iqr_time), axis=0)
    result_df = pd.concat([b, df.iloc[:, metab_num:]], axis=1)
    nums = result_df.iloc[:,:metab_num-len(abnormal_metabolites)].apply(lambda x : count_outliers(x,iqr_time), axis=0)
    # print(nums)
    fig.add_subplot(122)
    plt.hist(nums)
    return result_df
