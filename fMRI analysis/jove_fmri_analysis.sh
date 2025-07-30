#!/bin/bash
#
# Rename the input files to match:
# - in ${subj}.recon subdirectory:
#  epitmt1.nii.gz epitmt2.nii.gz t1anat.nii.gz  (fMRI data)
# - in ${subj}.physio subdirectory:
#  epitmt1.puls.1D epitmt1.resp.1D              (pulse and resp data)
#  epitmt2.puls.1D epitmt2.resp.1D
# - in ${subj}.paradigm subdirectory:
#  tmt-a.txt tmt-b.txt                          (stimulus timing [fixed])
#  tmt-a${subj}epirt.txt                        (completion times [measured])
#  tmt-b${subj}epirt.txt
#
# Requires: AFNI

##################################################

# How many timepoints to ignore and censor from beginning of each run
ignore=4
# Motion threshold for censoring
mthr=0.5

##################################################

# Subject ID
subj=$1
shift

#---------- Prepare stim files ----------

# Merge onset time with completion time for duration modulated analysis
(
cd $subj.paradigm

1dMarry tmt-a.txt tmt-a${subj}epirt.txt >tmt-a${subj}epidm.txt
1dMarry tmt-b.txt tmt-b${subj}epirt.txt >tmt-b${subj}epidm.txt
)

#---------- Prepare the physio data ----------

# Convert physio data into ricor regressors
(
cd $subj.physio

RetroTS.py -r epitmt1.resp.1D -c epitmt1.puls.1D -p 400 \
 -n `3dinfo -nk ../$subj.recon/epitmt1.nii.gz` \
 -v `3dinfo -tr ../$subj.recon/epitmt1.nii.gz` \
 -slice_order custom \
 -slice_offset "`3dinfo -slice_timing -sb_delim ' ' ../$subj.recon/epitmt1.nii.gz`" \
 -prefix epitmt1
RetroTS.py -r epitmt2.resp.1D -c epitmt2.puls.1D -p 400 \
 -n `3dinfo -nk ../$subj.recon/epitmt2.nii.gz` \
 -v `3dinfo -tr ../$subj.recon/epitmt2.nii.gz` \
 -slice_order custom \
 -slice_offset "`3dinfo -slice_timing -sb_delim ' ' ../$subj.recon/epitmt2.nii.gz`" \
 -prefix epitmt2
)

#---------- Crop anatomy ----------

# Zero out neck in the volume
(
cd $subj.recon

3dcalc -prefix t1anat_crop.nii.gz -a t1anat.nii.gz -expr 'ispositive(z+80)*a' -RAI
)

#---------- Anatomical skull stripping ----------

# Just using the final skull stripped output in orig space
(
cd $subj.recon

@SSwarper -input t1anat_crop.nii.gz -base MNI152_2009_template_SSW.nii.gz -subid ssw -SSopt '-blur_fwhm 3'
3dcopy anatSS.ssw.nii t1anat_crop_brain.nii.gz
)

#---------- Do the bulk of the work ----------

afni_proc.py -subj_id $subj                                            \
-out_dir $subj.results                                                 \
-script proc_subj_${subj}.sh -scr_overwrite -execute                   \
-blocks despike ricor tshift align tlrc volreg blur mask scale regress \
-dsets $subj.recon/epitmt[1-2].nii.gz                                  \
-copy_anat $subj.recon/t1anat_crop_brain.nii.gz                        \
-anat_has_skull no                                                     \
-ricor_regress_method across-runs                                      \
-ricor_regs $subj.physio/epitmt[1-2].slibase.1D                        \
-tlrc_base MNI152_2009_template_SSW.nii.gz                             \
-tlrc_opts_at -init_xform AUTO_CENTER                                  \
-tlrc_rmode quintic                                                    \
-tlrc_NL_warp                                                          \
-tlrc_NL_warped_dsets                                                  \
  $subj.recon/anatQQ.ssw.nii                                           \
  $subj.recon/anatQQ.ssw.aff12.1D                                      \
  $subj.recon/anatQQ.ssw_WARP.nii                                      \
-volreg_base_ind 1 $ignore                                             \
-volreg_align_e2a                                                      \
-volreg_opts_vr -twopass                                               \
-volreg_tlrc_warp                                                      \
-volreg_zpad 4                                                         \
-align_opts_aea -cost lpc+ZZ -giant_move -AddEdge -check_flip          \
-blur_size 5.0                                                         \
-blur_in_mask yes                                                      \
-mask_epi_anat yes                                                     \
-regress_censor_extern $subj.paradigm/epitmt-censorfirst7s+trial1.1D   \
-regress_censor_first_trs $ignore                                      \
-regress_censor_motion $mthr                                           \
-regress_censor_outliers 0.1                                           \
-regress_apply_mot_types demean deriv                                  \
-regress_motion_per_run                                                \
-regress_apply_ricor yes                                               \
-regress_stim_types AM1 AM1                                            \
-regress_stim_times $subj.paradigm/tmt-a${subj}epidm.txt               \
  $subj.paradigm/tmt-b${subj}epidm.txt                                 \
-regress_stim_labels tmta tmtb                                         \
-regress_basis_multi 'dmUBLOCK(1)' 'dmUBLOCK(1)'                       \
-regress_opts_3dD -jobs 4                                              \
 -gltsym 'SYM: -tmta tmtb'                                             \
 -glt_label 1 tmtb-tmta                                                \
 -gltsym 'SYM: 0.5*tmta 0.5*tmtb'                                      \
 -glt_label 2 avgtmt                                                   \
-mask_apply anat                                                       \
-regress_fout no                                                       \
-regress_no_iresp                                                      \
-regress_compute_fitts                                                 \
-regress_est_blur_errts                                                \
-regress_reml_exec

#---------- Cleanup ----------

# Save disk space
find $subj.results -type f -name \*.BRIK -exec gzip "{}"  \;

