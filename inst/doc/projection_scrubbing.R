## ----include=FALSE------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(autodep = TRUE, cache = FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
# devtools::install_github('mandymejia/fMRIscrub')
library(fMRIscrub)

## -----------------------------------------------------------------------------
dim(Dat1) # subject 0050048 (Pitt); many artifacts
dim(Dat2) # subject 0051479 (Caltech); few or no artifacts

## -----------------------------------------------------------------------------
ps.Dat1 = pscrub(Dat1, verbose=TRUE, comps_mean_dt=1, comps_var_dt=1)
# same as: scrub(Dat1, verbose=TRUE, comps_mean_dt=1, comps_var_dt=1)
ps.Dat2 = pscrub(Dat2, verbose=TRUE, comps_mean_dt=1, comps_var_dt=1)

## ----fig.width=8, fig.height=3------------------------------------------------
p1 <- plot(ps.Dat1, title="Dat1", show.legend=FALSE)
p2 <- plot(ps.Dat2, title="Dat2", show.legend=FALSE)
cowplot::plot_grid(p1, p2, nrow=1)

## ----fig.width=8, fig.height=3------------------------------------------------
p1 <- plot(DVARS(Dat1), title="Dat1", show.legend=FALSE)
# same as: plot(scrub(Dat1, "DVARS"), title="Dat1", show.legend=FALSE)
p2 <- plot(DVARS(Dat2), title="Dat2", show.legend=FALSE)
cowplot::plot_grid(p1, p2, nrow=1)

## ----warning=FALSE------------------------------------------------------------
library(oro.nifti)
library(neurobase)

## -----------------------------------------------------------------------------
# Unmask the first scan
fname = system.file("extdata", "Dat1_mask.nii.gz", package = "fMRIscrub")
Mask1 = readNIfTI(fname) > 0 #Pitt_0050048 (full of artifacts)
Mask1 = array(Mask1, dim=c(dim(Mask1), 1)) # 2D --> 3D slice
Img1 = fMRIscrub::unmask_vol(t(Dat1), Mask1)

## ---- fig.width=7, fig.height=4-----------------------------------------------
mfrow_original <- par("mfrow")
par(mfrow=c(1,2))
levs = ps.Dat1$measure
t_med = order(levs)[ceiling(length(levs)/2)]
t_max = which.max(levs)

image(Img1[,,,t_med], main=paste0('Median lev (T = ', t_med, ')'))
image(Img1[,,,t_max], main=paste0('Maximum lev (T = ', t_max, ')'))

## -----------------------------------------------------------------------------
par(mfrow=mfrow_original)

## ---- fig.width=7, fig.height=4-----------------------------------------------
psx <- pscrub(Dat1, projection="ICA", get_dirs=TRUE, comps_mean_dt=1, comps_var_dt=1)
artImg1 = artifact_images(psx)

par(mfrow=c(1,2))

# Constant voxels are deleted during the `pscrub` algorithm, so the artifact images will have
# missing values where the constant voxels were. 
artImg1.mean = unmask_vol(t(artImg1$mean), Mask1)
artImg1.top = unmask_vol(t(artImg1$top), Mask1)

idx = which(which(psx$outlier_flag) == t_max)
image(artImg1.mean[,,1,idx], main=paste0('Lev image, mean (T=',t_max,')'))
image(artImg1.top[,,1,idx], main=paste0('Lev image, top (T=',t_max,')'))

## ----fig.height=3, fig.width=6------------------------------------------------
# Try all projections: this is excessive to plot
# ps.Dat1All <- fMRIscrub:::pscrub_multi(
#   Dat1, projection="all", verbose=TRUE, comps_mean_dt=1, comps_var_dt=1
# )

#  the default (ICA + kurtosis), fusedPCA + kurtosis, and ICA
ps.Dat1.3 <- fMRIscrub:::pscrub_multi(
  Dat1, projection=c("ICA_kurt", "fusedPCA_kurt", "ICA"), verbose=TRUE, comps_mean_dt=1, comps_var_dt=1
)
fMRIscrub:::plot.scrub_projection_multi(ps.Dat1.3, legend.position="bottom")

