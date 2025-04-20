
import numpy as np
import pandas as pd

import plotly
import plotly.graph_objs
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from scipy import stats
from scipy.stats import ks_2samp

import matplotlib.pyplot as plt

BINS = 1000
z_score_max = 2
# px.colors.DEFAULT_PLOTLY_COLORS

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



def ks_test (wt, ko, x_label, xaxis_dtick): 
    
    D_auto, p_value = ks_2samp(wt, ko)
    
    all_values = np.concatenate((wt, ko))
    min_all_values = min(all_values)
    max_all_values = max(all_values)
    
    count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    
    count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
    pdf2 = count2 / sum(count2)
    cdf2 = np.cumsum(pdf2)
        
    D = abs(cdf-cdf2)
    maximum = max(D)
    
    x_maximum = bins_count[np.argmax(D)]
    
    y_min = cdf[np.argmax(D)]
    y_max = cdf2[np.argmax(D)]
    

    d = {'WT': cdf, 'KO': cdf2}
    df = pd.DataFrame(data=d)
    
    fig = px.line(df, x=bins_count[1:], y=['WT', 'KO'], template="simple_white", width=800, height=800, line_shape="vh")

    fig.add_trace(
        
        go.Scatter(
        x=[x_maximum, x_maximum],
        y=[y_min, y_max], 
        mode='lines', 
        line=dict(color='Black', dash='dash'), 
        showlegend=False))
    
    fig.update_layout(
        
        autosize= False, width=1000, height=1000, margin = dict(l=100, r=100, b=100, t=100, pad=4),    

        xaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = xaxis_dtick),         
        yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 0.2),  
        
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.2), 
        legend_title_text=None, 
        
        #title="Plot Title",
        xaxis_title=x_label,
        yaxis_title="CDF",

        font=dict(family="arial", size=36, color="Black"))
    
    plotly.offline.plot(fig, filename= "html_kstest_" +x_label + ".html")
    plt.figure() 
    
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

D_area_xx, D_area_xx_manual_calc, p_value_area_xx = ks_test (area_xx_plus, area_xx_minus, "Relative territory Chr. #xx", 0.1)
D_area_yy, D_area_yy_manual_calc, p_value_area_yy = ks_test (area_yy_plus, area_yy_minus, "Relative territory Chr. #yy", 0.1)
D_circ_xx, D_circ_xx_manual_calc, p_value_circ_xx = ks_test (circ_xx_plus, circ_xx_minus, "Circularity of Chr. #xx", 0.5)
D_circ_yy, D_circ_yy_manual_calc, p_value_circ_yy = ks_test (circ_yy_plus, circ_yy_minus, "Circularity of Chr. #yy", 0.5)
D_coloc_xx_yy, D_coloc_xx_yy_manual_calc, p_value_coloc_xx_yy = ks_test (coloc_xx_yy_plus, coloc_xx_yy_minus, "Relative colocalization", 0.05)






def ks_test_subplot (all_data, all_xlabels): 
    
    fig = make_subplots(rows=1, cols=5, shared_yaxes=True) #,horizontal_spacing=0.05, vertical_spacing=0.1

    for i in range(5):

        wt = all_data[i][0]
        ko = all_data[i][1]

        D_auto, p_value = ks_2samp(wt, ko)
        
        all_values = np.concatenate((wt, ko))
        min_all_values = min(all_values)
        max_all_values = max(all_values)
        
        count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        
        count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
        pdf2 = count2 / sum(count2)
        cdf2 = np.cumsum(pdf2)

        D = abs(cdf-cdf2)
        maximum = max(D)
        
        x_maximum = bins_count[np.argmax(D)]
        
        y_min = cdf[np.argmax(D)]
        y_max = cdf2[np.argmax(D)]
  
        fig.add_trace(
            
            go.Scatter(
            x= bins_count[1:],
            y= cdf, 
            name= 'WT',
            line= dict(color='rgb(31, xx9, 180)'),
            showlegend= i==0), 
            row=1, col=i+1)            
        
        fig.add_trace(
            
            go.Scatter(            
            x= bins_count[1:],
            y= cdf2, 
            name= 'KO',
            line= dict(color='rgb(255, 127, yy)'),
            showlegend= i==0), 
            row=1, col=i+1)           

        
        fig.add_trace(
            
            go.Scatter(
            x= [x_maximum, x_maximum],
            y= [y_min, y_max], 
            mode= 'lines', 
            line= dict(color='Black', dash='dash'), 
            showlegend= False), 
            row=1, col=i+1)            
        
        if i ==0:
            fig['layout']['xaxis']['title'] = all_xlabels[i]
        else:
            fig['layout']['xaxis'+str(i+1)]['title'] = all_xlabels[i]
        
    fig.update_layout(
        
        autosize=False, width=2000, height=565, margin = dict(l=100, r=100, b=100, t=100, pad=4),    
        
        xaxis= dict(mirror=True, ticks='outside', showline=True),        
        yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 0.5),         

        legend_title_text=None,
        #legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.2),
         
        #title="Plot Title",
        #xaxis_title="x_title",
        yaxis_title="CDF",
        
        font=dict(family="arial", size=36, color="Black"),
        
        template="simple_white")
    
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)

    plotly.offline.plot(fig, filename= "html_kstest_subplot.html")


all_data = [(area_xx_plus, area_xx_minus),
            (area_yy_plus, area_yy_minus),
            (circ_xx_plus, circ_xx_minus),
            (circ_yy_plus, circ_yy_minus),
            (coloc_xx_yy_plus, coloc_xx_yy_minus)] 
             
    
all_xlabels = ["Relative territory<br>Chr. #xx", 
               "Relative territory<br>Chr. #yy", 
               "Circularity of<br>Chr. #xx", 
               "Circularity of<br>Chr. #yy", 
               "Relative<br>colocalization"]


ks_test_subplot (all_data, all_xlabels)

