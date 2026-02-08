import anndata
import scanpy as sc
import scipy as s
from scipy.sparse import csr_matrix, issparse

def load_adata(adata_file, metadata_file = None, normalise = False, cells = None, cell_column = "cell", features = None, filter_lowly_expressed_genes = False, set_colors = False, keep_counts=False):

	adata = sc.read(adata_file)

	# Convert to sparse matrices
	if not s.sparse.issparse(adata.X):
		adata.X = csr_matrix(adata.X)
	if len(adata.layers.keys())>0:
		for i in list(adata.layers.keys()):
			if not issparse(adata.layers[i]):
				adata.layers[i] = csr_matrix(adata.layers[i])

	if cells is not None:
		tmp = np.mean(np.isin(cells,adata.obs.index.values)==False)
		if tmp<1: print("%.2f%% of cells provided are not observed in the adata, taking the intersect..." % (100*tmp))
		cells = np.intersect1d(cells,adata.obs.index.values)
		adata = adata[cells,:]

	if features is not None:
		adata = adata[:,features]

	if metadata_file is not None:
		metadata = pd.read_table(metadata_file, delimiter="\t", header=0).set_index(cell_column, drop=False)
		metadata = metadata.loc[cells]
		assert np.all(adata.obs.index.isin(metadata[cell_column]))
		# assert np.all(metadata.cell.isin(adata.obs.index))
		assert metadata.shape[0] == adata.shape[0]
		adata.obs = metadata#.reindex(adata.obs.index)

	if filter_lowly_expressed_genes:
		sc.pp.filter_genes(adata, min_counts=10)

	if keep_counts:
		adata.layers["raw"] = adata.X.copy()

	if normalise:
		sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=False)
		sc.pp.log1p(adata)

	if set_colors:
		colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['celltype']))]
		adata.uns['celltype_colors'] = colPalette_celltypes
		colPalette_stages = [opts["stage_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['stage']))]
		adata.uns['stage_colors'] = colPalette_stages

	return adata

def scale(X, x_min, x_max):
    nom = (X - X.min(axis=0)) * (x_max - x_min)
    denom = X.max(axis=0) - X.min(axis=0)
    denom[denom == 0] = 1
    return x_min + nom / denom


# cmap = custom_div_cmap(11, mincol='g', midcol='0.9' ,maxcol='CornflowerBlue')
def custom_div_cmap(numcolors=11, name='custom_div_cmap',
                    mincol='blue', midcol='white', maxcol='red'):
    """ 
    Default is blue to white to red with 11 colors.  
    Colors can be specified in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap 
    cmap = LinearSegmentedColormap.from_list(name=name, colors =[mincol, midcol, maxcol], N=numcolors)
    return cmap