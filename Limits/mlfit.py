#!/usr/bin/env python
# coding: utf-8

# In[110]:


from ROOT import TFile
from ROOT import gDirectory
from ROOT import TCanvas
from ROOT import TGraphErrors
from ROOT import TMultiGraph
from ROOT import TLegend


# In[111]:


mlfit = TFile('mlfit.root', 'READ')
mlfit.ls()


# In[149]:


prefit_hists = {'total_background':[], 'ewk':[], 'qcd':[], 'zgg':[]}
postfit_hists = {'total_background':[], 'ewk':[], 'qcd':[], 'zgg':[]}

for i in [1,2,3,4,5,6]:

    mlfit.cd('shapes_prefit/ch1_bin'+str(i))

    prefit_hists['total_background'].append(gDirectory.Get('total_background;1'))
    prefit_hists['ewk'].append(gDirectory.Get('ewk;1'))
    prefit_hists['qcd'].append(gDirectory.Get('qcd;1'))
    prefit_hists['zgg'].append(gDirectory.Get('zgg;1'))
    
    mlfit.cd('shapes_fit_b/ch1_bin'+str(i))
    
    postfit_hists['total_background'].append(gDirectory.Get('total_background;1'))
    postfit_hists['ewk'].append(gDirectory.Get('ewk;1'))
    postfit_hists['qcd'].append(gDirectory.Get('qcd;1'))
    postfit_hists['zgg'].append(gDirectory.Get('zgg;1'))


# In[150]:


from array import array

channel = 'total_background'
#channel = 'ewk'
#channel = 'qcd'
#channel = 'zgg'

X = array('f')
Xerr = array('f')

Ypre = array('f')
Ypre_err = array('f')

Ypost = array('f')
Ypost_err = array('f')

for i in range(0,6):
    
    X.append(i+1)
    Xerr.append(0.5)
    
    hpre = prefit_hists[channel][i]
    
    Ypre.append(hpre.GetBinContent(1))
    Ypre_err.append(hpre.GetBinError(1))

    hpost = postfit_hists[channel][i]
    
    Ypost.append(hpost.GetBinContent(1))
    Ypost_err.append(hpost.GetBinError(1))
    
    
c1 = TCanvas()    

gr_pre = TGraphErrors(6, X, Ypre, Xerr, Ypre_err)
gr_pre.SetMarkerColor(4)
gr_pre.SetMarkerStyle(21)
gr_pre.SetLineColor(4)

gr_post = TGraphErrors(6, X, Ypost, Xerr, Ypost_err)
gr_post.SetMarkerColor(2)
gr_post.SetMarkerStyle(21)
gr_post.SetLineColor(2)

mg = TMultiGraph()
mg.Add(gr_pre, 'p')
mg.Add(gr_post, 'p')

mg.SetTitle(channel)
mg.SetTitle(channel+';bin')

leg = TLegend(0.7, 0.7, 0.9, 0.9) 
leg.SetFillColor(0) 
leg.AddEntry(gr_pre, 'prefit', 'p')
leg.AddEntry(gr_post, 'postfit', 'p')

mg.Draw('a')
leg.Draw()

c1.Draw()


# In[151]:


mlfit.cd('shapes_prefit')
total_pre = gDirectory.Get('total_background')
total_pre.SetLineColor(4)

mlfit.cd('shapes_fit_b')
total_post = gDirectory.Get('total_background')
total_post.SetLineColor(2)

c2 = TCanvas()    

total_pre.Draw()
total_post.Draw('same')

leg = TLegend(0.7, 0.7, 0.9, 0.9) 
leg.SetFillColor(0) 
leg.AddEntry(total_pre, 'prefit', 'l')
leg.AddEntry(total_post, 'postfit', 'l')

leg.Draw()

c2.Draw()


# In[ ]:




