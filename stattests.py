from numpy import *
import matplotlib
import pylab
from numpy.random import randn
import numpy as numpy
import sys
from cluster import cluster
from scipy.stats import kstest
from scipy.stats import anderson

## Threshold defined by gmeans
ad_thresh = 1.0

ks_thresh = 0.05

dist = linspace(0,4,10)
ADstat = []
KSstat = []
ADmean = []
KSmean = []
ADhist = []
KShist = []
## Plots the normalized KS and AD test statistics as a function of separation distance
## between two distributions
trials = 100
members = 100
for i in dist:
    x = []
    y = []
    for j in range(trials):
        x1 = randn(members)
        y1 = randn(members)

        x2 = randn(members)+i
        y2 = randn(members)

        x = append(x1,x2)
        y = append(y1,y2)

        data_list = []

        for m in range(len(x)):
            data = array([x[m],y[m]])
            data_list.append(data)

        data_list = vstack(data_list)
        ks = kstest(x,'norm')
        ks1 = kstest(y1,'norm')
        ks,ad,cent,labels = cluster(data_list)
        KSstat = append(KSstat,ks)
        ADstat = append(ADstat,ad)
    
    AD_mu = mean(ADstat)
    KS_mu = mean(KSstat)
    
    ADmean = append(ADmean,AD_mu)
    KSmean = append(KSmean,KS_mu)
    
ADmean = ADmean/max(abs(ADmean))
KSmean = KSmean/max(abs(KSmean))

pylab.figure()
pylab.subplot(2,1,1)
pylab.plot(dist,ADmean,label='AD')
pylab.subplot(2,1,2)
pylab.plot(dist,KSmean,label='KS')
pylab.legend()
pylab.xlabel('Distance from Centers')
pylab.ylabel('Kolgomorov-Smirnov Test Statistic')
pylab.savefig('KSADtestvalues.png')
pylab.close()



## Plots the contrasted heat maps for the KS pvalue and the AD test statistic
## as the separation distance and cluster ratios vary
res = 100
ref_members = 1000
dist = linspace(-5,5,res)
gridks = zeros((res,len(dist)))
gridad = zeros((res,len(dist)))
m = 0
n = 0
ratio = linspace(ref_members,5*ref_members,res)
for k in dist:
    for r in ratio:
        x1 = randn(ref_members)-k/2
        x2 = randn(r)+k/2
        x = append(x1,x2)
        x = (x-mean(x))/std(x)
        x = sort(x)
        ks = kstest(x,'norm')
        ad = anderson(x,'norm')
        if ks[1] > ks_thresh:
            gridks[n,m] = 1
        elif ks[1] <= ks_thresh:
            gridks[n,m] = 0
        if ad[0] < ad_thresh:
            gridad[n,m] = 1
        elif ad[0] >= ad_thresh:
            gridad[n,m] = 0
        n = n + 1
    
    n = 0
    m = m + 1

pylab.subplot(1,2,1)
pylab.imshow((gridks),origin='lower')
pylab.title('KS Test')
pylab.xlabel('Separation Distance')
pylab.ylabel('Cluster Ratio')
pylab.xticks(linspace(1,100,5),linspace(-5,5,5))
pylab.yticks(linspace(1,100,5),arange(5)+1)
pylab.subplot(1,2,2)
pylab.imshow((gridad),origin='lower')
pylab.title('AD Test')
pylab.xlabel('Separation Distance')
pylab.ylabel('Cluster Ratio')
pylab.xticks(linspace(1,100,5),linspace(-5,5,5))
pylab.yticks(linspace(1,100,5),arange(5)+1)
pylab.savefig('heat_map.png')
pylab.close()

## Checks distributions of 2-dimensional clusters calculated KS pvalues and AD test statistic values
## by using the cluster.py script while separation distance can be varied
ad2m = []
ks2m = []
ad1m = []
ks1m = []
fpks = []
fpad = []
fnks = []
fnad = []
stdev = 1.5
offset = 100
members = 50
trials = 10000


for t in range(trials):
    x = []
    ## initialize single distribution
    x0 = randn(2*j)+offset
    x0 = x0/max(x0)
    x0 = (x0-mean(x0))/std(x0)
    ## initialize double distributions
    x1 = randn(j)
    sep = stdev*numpy.std(x1)
    x2 = randn(j) + stdev*sep
    x = append(x1,x2)+offset
    x = x/max(x)
    x = (x-mean(x))/std(x)

    ## kstest
    ks1d = kstest(x0,'norm')
    ks2d = kstest(x,'norm')
    ks1m = append(ks1m,ks1d[1])
    ks2m = append(ks2m,ks2d[1])
    ## adtest
    ad1d = anderson(x0,'norm')
    ad2d = anderson(x,'norm')
    ad1m = append(ad1m,ad1d[0])
    ad2m = append(ad2m,ad2d[0])


#test = pylab.hist(ksm,bins=50)
#print test
pylab.figure()
pylab.title('KS and AD Value Distribution Comparisons')
pylab.subplot(2,1,1)
pylab.hist(ad1m,bins=50,label='Single Distribution')
pylab.hist(ad2m,bins=50,label='2 Distributions %i-STD apart' %round(stdev))
pylab.legend()
pylab.title('AD Statistic Gaussian/NonGaussian Comparison')
pylab.xlabel('AD Statistic')
pylab.ylabel('Frequency of AD Statistic')

pylab.subplot(2,1,2)
pylab.hist(10*ks1m,bins=50,label = 'Single Distribution')
pylab.hist(10*ks2m,bins=50,label='2 Distributions %i-STD apart' %round(stdev))
pylab.title('KS Statistic Gaussian/NonGaussian Comparison')
pylab.xlabel('KS Statistic')
pylab.ylabel('Frequency of KS Statistic')
pylab.legend()
pylab.savefig('ksaddistcompare.png')
pylab.show()
pylab.close()

print 'Anderson Darling'
tad = ad1m>ad_thresh
pacvad = 100*len(ad1m[tad == True])/float(len(ad1m))
print 'False Negatives(AD)',pacvad
fad = ad2m<ad_thresh
pbcvad = 100*len(ad2m[fad == True])/float(len(ad2m))
print 'False Positives(AD)',pbcvad
print 'Kolmogorov-Smirnov'
fks = ks1m<ks_thresh
pacvks = 100*len(ks1m[fks == True])/float(len(ks1m))
print 'False Negatives(KS)',pacvks
tks = ks2m>ks_thresh
pbcvks = 100*len(ks2m[tks == True])/float(len(ks2m))
print 'False Positives(KS)',pbcvks

#ad_single = pylab.hist(ad1m,bins=50)
#ad_double = pylab.hist(adm,bins=50)
#for i in range(50):
#    for j in range(50):
#        ads = ad_single[1]
#        add = ad_double[1]
#        if ad_single[0][i] == ad_double[0][j]:
#            if round(ads[i],3) == round(add[j],3):
#                print ads[i]

#ks_single = pylab.hist(ks1m,bins=50)
#ks_double = pylab.hist(ksm,bins=50)
#for i in range(50):
#    for j in range(50):
#        kss = ks_single[1]
#        ksd = ks_double[1]
#        if ks_single[0][i] == ks_double[0][j]:
#            if round(kss[i],3) == round(ksd[j],3):
#                print kss[i]


