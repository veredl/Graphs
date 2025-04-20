
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
from scipy.stats import ks_2samp

sns.set(rc={"figure.dpi":300, 'savefig.dpi':300, "figure.figsize":(4,4)})
sns.set_context('notebook')
sns.set_style("ticks")

BINS = 1000
z_score_max = 2


"""--------------------area--------------------"""

area_xx_minus = pd.read_csv('area_xx_minus.txt')
area_xx_plus = pd.read_csv('area_xx_plus.txt')
area_yy_minus = pd.read_csv('area_yy_minus.txt')
area_yy_plus = pd.read_csv('area_yy_plus.txt')

area_xx_minus = area_xx_minus[(np.abs(stats.zscore(area_xx_minus)) < z_score_max).all(axis=1)]
area_xx_plus = area_xx_plus[(np.abs(stats.zscore(area_xx_plus)) < z_score_max).all(axis=1)]
area_yy_minus = area_yy_minus[(np.abs(stats.zscore(area_yy_minus)) < z_score_max).all(axis=1)]
area_yy_plus = area_yy_plus[(np.abs(stats.zscore(area_yy_plus)) < z_score_max).all(axis=1)]


"""--------------------circ--------------------"""

circ_xx_minus = pd.read_csv('circ_xx_minus.txt')
circ_xx_plus = pd.read_csv('circ_xx_plus.txt')
circ_yy_minus = pd.read_csv('circ_yy_minus.txt')
circ_yy_plus = pd.read_csv('circ_yy_plus.txt')

circ_xx_minus = circ_xx_minus[(np.abs(stats.zscore(circ_xx_minus)) < z_score_max).all(axis=1)]
circ_xx_plus = circ_xx_plus[(np.abs(stats.zscore(circ_xx_plus)) < z_score_max).all(axis=1)]
circ_yy_minus = circ_yy_minus[(np.abs(stats.zscore(circ_yy_minus)) < z_score_max).all(axis=1)]
circ_yy_plus = circ_yy_plus[(np.abs(stats.zscore(circ_yy_plus)) < z_score_max).all(axis=1)]


"""--------------------coloc--------------------"""

coloc_xx_yy_minus = pd.read_csv('coloc_xx_yy_minus.txt')
coloc_xx_yy_plus = pd.read_csv('coloc_xx_yy_plus.txt')

coloc_xx_yy_minus = coloc_xx_yy_minus[(np.abs(stats.zscore(coloc_xx_yy_minus)) < z_score_max).all(axis=1)]
coloc_xx_yy_plus = coloc_xx_yy_plus[(np.abs(stats.zscore(coloc_xx_yy_plus)) < z_score_max).all(axis=1)]



def ks_test (wt, ko, x_label): 
    
    D_auto, p_value = ks_2samp(wt, ko)
    
    all_values = np.concatenate((wt, ko))
    min_all_values = min(all_values)
    max_all_values = max(all_values)
    
    count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    plt.step(bins_count[1:], cdf, label="WT")
    
    count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
    pdf2 = count2 / sum(count2)
    cdf2 = np.cumsum(pdf2)
    plt.step(bins_count2[1:], cdf2, label="KO") # bins_count2 is the same as bins_count
        
    D = abs(cdf-cdf2)
    maximum = max(D)
    
    x_maximum = bins_count[np.argmax(D)]
    
    y_min = cdf[np.argmax(D)]
    y_max = cdf2[np.argmax(D)]
    
    plt.plot((x_maximum, x_maximum), (y_min, y_max), 'k' , linestyle="--")
    plt.xlabel(x_label)
    plt.ylabel("CDF")

    plt.legend(frameon=False)

    plt.show()
    
    return D_auto, maximum , p_value


area_xx_plus = area_xx_plus.to_numpy().flatten()
area_xx_minus = area_xx_minus.to_numpy().flatten()

area_yy_plus = area_yy_plus.to_numpy().flatten()
area_yy_minus = area_yy_minus.to_numpy().flatten()

circ_xx_plus = circ_xx_plus.to_numpy().flatten()
circ_xx_minus = circ_xx_minus.to_numpy().flatten()

circ_yy_plus = circ_yy_plus.to_numpy().flatten()
circ_yy_minus = circ_yy_minus.to_numpy().flatten()

coloc_xx_yy_plus = coloc_xx_yy_plus.to_numpy().flatten()
coloc_xx_yy_minus = coloc_xx_yy_minus.to_numpy().flatten()

D_area_xx, D_area_xx_manual_calc, p_value_area_xx = ks_test (area_xx_plus, area_xx_minus, "Relative territory Chr. #xx")
D_area_yy, D_area_yy_manual_calc, p_value_area_yy = ks_test (area_yy_plus, area_yy_minus, "Relative territory Chr. #yy")
D_circ_xx, D_circ_xx_manual_calc, p_value_circ_xx = ks_test (circ_xx_plus, circ_xx_minus, "Circularity of Chr. #xx")
D_circ_yy, D_circ_yy_manual_calc, p_value_circ_yy = ks_test (circ_yy_plus, circ_yy_minus, "Circularity of Chr. #yy")
D_coloc_xx_yy, D_coloc_xx_yy_manual_calc, p_value_coloc_xx_yy = ks_test (coloc_xx_yy_plus, coloc_xx_yy_minus, "Relative colocalization")



def ks_test_subplot (all_data, all_xlabels): 
    
    fig, axs = plt.subplots(nrows = 1, ncols = 5, sharey=True, figsize=(9, 2))
    fig.tight_layout()    
    sns.set_theme(style="ticks",font="arial") #,fontsize=28)
   
    for i in range(5):
          
        wt = all_data[i][0]
        ko = all_data[i][1]
        x_label = all_xlabels[i]
        
        D_auto, p_value = ks_2samp(wt, ko)
        
        all_values = np.concatenate((wt, ko))
        min_all_values = min(all_values)
        max_all_values = max(all_values)
        
        count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        axs[i].plot(bins_count[1:], cdf, label="WT")
        
        count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
        pdf2 = count2 / sum(count2)
        cdf2 = np.cumsum(pdf2)
        axs[i].plot(bins_count2[1:], cdf2, label="KO") # bins_count2 is the same as bins_count  
        
        D = abs(cdf-cdf2)
        maximum = max(D)
        
        x_maximum = bins_count[np.argmax(D)]
        
        y_min = cdf[np.argmax(D)]
        y_max = cdf2[np.argmax(D)]
            
        axs[i].plot((x_maximum, x_maximum), (y_min, y_max), 'k' , linestyle="--")     
        axs[i].set_xlabel(x_label)
   
        if i == 0:
            axs[i].set_ylabel('CDF')
    
        plt.legend(frameon=False)
    

sns.set()
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_context('notebook')
sns.set_style("ticks")


all_data = [(area_xx_plus, area_xx_minus),
            (area_yy_plus, area_yy_minus),
            (circ_xx_plus, circ_xx_minus),
            (circ_yy_plus, circ_yy_minus),
            (coloc_xx_yy_plus, coloc_xx_yy_minus)] 
             
    
all_xlabels = ["Relative territory\nChr. #xx", 
               "Relative territory\nChr. #yy", 
               "Circularity of\nChr. #xx", 
               "Circularity of\nChr. #yy", 
               "Relative\ncolocalization"]


ks_test_subplot (all_data, all_xlabels)
