import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import gzip
import numpy as np
import pandas as pd
import pymc3 as pm
import theano.tensor as tt
import toyplot
from statsmodels.stats.moment_helpers import corr2cov, cov2corr
import h5py
#import dask.array as da #maybe not necessary?



# opening a file handle
with h5py.File("db.h5", 'w') as io5:
    
    # initializing an empty h5 database
    nrows = 24006
    ncols = 5903
    
    db = io5.create_dataset(name='scCounts', shape=(nrows, ncols), dtype=np.float64, compression='gzip')

    # iterate over chunks of the input data and fill in the db
    ## grab ~4K lines of data
    x = 6
    chunks = 3999
    
    bigdata = gzip.open("C:\\Users\\Alexander\\PDSB\\project\\GSE103322_HNSCC_all_data.txt.gz", 'r')
    for chunk in range(x, nrows, chunks):

        # fill the array
        listdata = []
        for line in range(chunk):
            listdata.append(bigdata.readline().decode().split("\t"))
        
        # convert list of 20K lines to array of floats
        data = np.array(listdata)#.astype(np.float64)
        
        # stick into db
        db[chunk:chunk+chunks, :] = data
        
    bigdata.close()


#saving gene names
names = pd.read_csv('C:\\Users\\Alexander\\PDSB\\project\\GSE103322_HNSCC_all_data.txt.gz',
                    sep="\t", usecols = [0], squeeze=True)



glis = ['SPARC', 'KLF4', 'MALAT1', #EMT, my data
        'CCR7', 'CXCR2',  #HNSCC specific biomarkers, apparently
        'MMP9',  #metastasis
        'VEGFA', #angiogenesis
        'KRAS', 'MYC',  #proliferation
        'MIF'] #immune penetration

#regex magic
searchforthese = ['\\b{0}'.format(i) for i in glis]
searchforthese = ['{0}\\b'.format(i) for i in searchforthese]
searchforthese

searchindex = []
for i in searchforthese:
    print(i)
    print(names[names.str.contains(i) == 1])
    print('\n')
    searchindex.append(names[names.str.contains(i) ==1].index[0])


#getting data
cols = range(1,5903,13)
gem = np.loadtxt('C:\\Users\\Alexander\\PDSB\\project\\GSE103322_HNSCC_all_data.txt.gz',
	skiprows=6, usecols = cols)    

x = gem[searchindex,:]
xsize = x.shape[1]
xmeans = x.mean(axis=0)


#using pymc3
with pm.Model() as model:
    
    # distribution to draw stds from (not a rv for now) but is a wide prior
    stds = pm.HalfCauchy.dist(2)
    
    # sample correlations by drawing halfCauchy std values
    # eta=1 is uniform corrs, eta>1 puts more weight on corr=0.
    packed = pm.LKJCholeskyCov('packed_chol', n=xsize, eta=1., sd_dist=stds)
    
    # expand correlation matrix and define VCV as dot product
    chol = pm.expand_packed_triangular(xsize, packed)
    cov = pm.Deterministic('cov', chol.dot(chol.T))
    sigma = pm.Deterministic('sigma', tt.sqrt(tt.diag(cov)))

    # draw the means from a normal with starting values from observed means
    # set this mean and std of this dist according to your data, or it could
    # be made into a random variable to be fit.
    means = pm.Normal('means', 5, 5, shape=xsize, testval=xmeans)
    obs = pm.MvNormal('obs', means, chol=chol, observed=x)

with model:
    trace = pm.sample(tune=1000, draws=5000, njobs=3, cores=3)