# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 22:25:24 2019

@author: Wang Aki
"""
import os
os.chdir("D:\\A_My_Data\\HCC\\1_data\\validation_2")

import numpy as np
#from sklearn import metrics
from scipy import interp
#import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
#from matplotlib.backends.backend_pdf import PdfPages
import scipy.io as scio
import pandas as pd 



try:
    import matplotlib.pyplot as plt
except:
    raise
import networkx as nx
#用grid_2d_graph()生成一个16个节点的网格图
"""
G=nx.grid_2d_graph(4,4)  #4x4 grid
pos=nx.spring_layout(G,iterations=100)
#開始画各个小图
plt.subplot(221)
nx.draw(G,pos,font_size=8)
plt.subplot(222)
nx.draw(G,pos,node_color='k',node_size=0,with_labels=False)
plt.subplot(223)
nx.draw(G,pos,node_color='g',node_size=250,with_labels=False,width=6)
#最后一幅子图转为有向图
plt.subplot(224)
H=G.to_directed()
nx.draw(H,pos,node_color='b',node_size=20,with_labels=False)
plt.savefig("four_grids.png")
plt.show()
"""




#dataFile = 'D:/shx_bioinformatics/TCGA_Project/used_for_plot_ROC_and_static/reg_s_33_new.mat'
#data = scio.loadmat(dataFile)
#Normal =data['Normal_pr']
#Normal.shape[1]
#Disease =data['D']
#d_used =Disease[:,1]#只有一个
##
##贴标签
#lable=np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
test_fitted = pd.read_table("v2_test_fitted.txt")
test_if_health = pd.read_table("v2_test_if_health.txt")

test_fitted_2 = pd.read_table("v2_test_fitted_2.txt")
test_if_health_2 = pd.read_table("v2_test_if_health_2.txt")

test_fitted_3 = pd.read_table("v2_test_fitted_3.txt")
test_if_health_3 = pd.read_table("v2_test_if_health_3.txt")




mean_fpr = np.linspace(0, 1, 1000)
tprs=[]
aucs=[]

#pdf = PdfPages('my_figure.pdf')
f=plt.figure(figsize=(20,25))


for i in range(0,20):
    mean_fpr = np.linspace(0, 1, 1000)
    tprs=[]
    aucs=[]
    plt.subplot(5,4,i+1)
    
    #probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test]) #模型训练，模型测试
    #.predict_proba：返回预测属于某标签的概率  #.predict：返回预测标签 
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(test_if_health.iloc[:,i], test_fitted.iloc[:,i])
    tprs.append(interp(mean_fpr, fpr, tpr)) #interp：求插值，fpr=mean_fpr时的tpr，添加在空list  tprs中
    tprs[-1][0] = 0.0  #第1个和最后一个元素赋值为0
    roc_auc = auc(fpr, tpr) #求auc值
    aucs.append(roc_auc) #auc值添加在空list  aucs中
    plt.plot(fpr, tpr, lw=1, alpha=0.3,  #画出roc曲线
             label='GSE64041 AUC = %0.2f' % ( roc_auc))
    fpr, tpr, thresholds = roc_curve(test_if_health_2.iloc[:,i], test_fitted_2.iloc[:,i])
    tprs.append(interp(mean_fpr, fpr, tpr)) #interp：求插值，fpr=mean_fpr时的tpr，添加在空list  tprs中
    tprs[-1][0] = 0.0  #第1个和最后一个元素赋值为0
    roc_auc = auc(fpr, tpr) #求auc值
    aucs.append(roc_auc) #auc值添加在空list  aucs中
    
    plt.plot(fpr, tpr, lw=1, alpha=0.3,  #画出roc曲线
             label='GSE45436 AUC = %0.2f' % ( roc_auc))
    
    fpr, tpr, thresholds = roc_curve(test_if_health_3.iloc[:,i], test_fitted_3.iloc[:,i])
    tprs.append(interp(mean_fpr, fpr, tpr)) #interp：求插值，fpr=mean_fpr时的tpr，添加在空list  tprs中
    tprs[-1][0] = 0.0  #第1个和最后一个元素赋值为0
    roc_auc = auc(fpr, tpr) #求auc值
    aucs.append(roc_auc) #auc值添加在空list  aucs中
    plt.plot(fpr, tpr, lw=1, alpha=0.3,  #画出roc曲线
             label='RNA_Seq AUC = %0.2f' % ( roc_auc))
    
    plt.plot([0, 1], [0, 1], linestyle='--', lw=3, color='darkorange',
          alpha=.8)#label='Chance',
    mean_tpr = np.mean(tprs, axis=0)#列求和
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='mediumpurple',
           label=r'Mean AUC = %0.2f $\pm$ %0.2f' % (mean_auc, std_auc),
           lw=3, alpha=.8)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2
                 ,label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('1-Specificity',size=13,x=0.5,y=0.35)
    plt.ylabel('Sensitivity',size=13,x=0.35,y=0.5)
    plt.title(r'M%d'%(i+1),size=15,y=1.0)
    
    plt.legend(loc="lower right")
    
    plt.tight_layout()# pad=0.2, h_pad=0.1, w_pad=0.5, rect=None)

    #plt.show()
#plt.show()    
 
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
#                wspace=0.5, hspace=1)

plt.show()   
#f.savefig("my_mean_roc.pdf", bbox_inches='tight')
f.savefig("v2_my_all_roc.pdf", bbox_inches='tight')
    


   
"""
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
         label='Chance', alpha=.8)
mean_tpr = np.mean(tprs, axis=0)#列求和
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)
std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('My_data_ROC_curve')
plt.legend(loc="lower right")
plt.show()
"""

#pdf.savefig()  
#plt.close()
#pdf.close()
