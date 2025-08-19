library(subtyper)
# =============================================================================
#
#   Part 2: The Comprehensive Test Suite
#
# =============================================================================

if (requireNamespace("testthat", quietly = TRUE)) {
  
  library(testthat)

  context("Testing ANTsPyMM Decoder - Structured Output")

  test_that("Handles edge cases and bad inputs", {
    res <- decode_antspymm_label(NULL)
    expect_equal(res$anatomy, "Invalid Input")
    res <- decode_antspymm_label("   ")
    expect_equal(res$anatomy, "Invalid Input")
  })
  
  test_that("Correctly decodes complex, squished tokens", {
    res_cit <- decode_antspymm_label("T1Hier_thk_mtg_sn_snc_leftcit168")
    expect_equal(res_cit$modality, "T1")
    expect_equal(res_cit$laterality, "Left")
    expect_equal(res_cit$anatomy, "Substantia Nigra Compacta")
    
    res_cereb <- decode_antspymm_label("T1Hier_vol_l_crus_icerebellum")
    expect_equal(res_cereb$modality, "T1")
    expect_equal(res_cereb$laterality, "Left")
    expect_equal(res_cereb$anatomy, "Cerebellum Crus I")
  })
  
  test_that("Correctly decodes DTI metrics", {
    res_fa <- decode_antspymm_label("DTI_mean_fa.anterior_corona_radiata.left.jhu_icbm_labels_1mm")
    expect_equal(res_fa$modality, "DTI")
    expect_equal(res_fa$laterality, "Left")
    expect_equal(res_fa$measurement, "Fractional Anisotropy")
    expect_equal(res_fa$anatomy, "Anterior Corona Radiata")
  })

  test_that("Correctly decodes rs-fMRI connectivity patterns", {
    res <- decode_antspymm_label("rsfMRI_fcnxpro122_DefaultA_2_Striatum")
    expect_equal(res$modality, "rs-fMRI")
    expect_equal(res$measurement, "Connectivity")
    expect_equal(res$anatomy, "Default Mode Network A to Striatal Network")
  })
  
  test_that("Correctly decodes QC and file metrics", {
    res_fn <- decode_antspymm_label("rsffn1")
    expect_equal(res_fn$modality, "rs-fMRI")
    expect_equal(res_fn$anatomy, "Filename")
    
    res_id <- decode_antspymm_label("flairid")
    expect_equal(res_id$modality, "T2-FLAIR")
    expect_equal(res_id$anatomy, "Image ID")
  })

  test_that("Correctly decodes common abbreviations and new patterns", {
    expect_equal(decode_antspymm_label("bn str cadp")$anatomy, "Caudate Nucleus")
    expect_equal(decode_antspymm_label("postcent_ctx")$anatomy, "Postcentral Cortex")
    expect_equal(decode_antspymm_label("postcent ctx")$anatomy, "Postcentral Cortex")
    expect_equal(decode_antspymm_label("T1w_mean_postcent_ctx")$anatomy, "Postcentral Cortex")
    expect_equal(decode_antspymm_label("T1Hier_vol_bn_str_cadp")$anatomy, "Caudate Nucleus")
    expect_equal(decode_antspymm_label("T1Hier_mean_sncdp")$anatomy, "Substantia Nigra Compacta")
    expect_equal(decode_antspymm_label("dti_fa_caudal_ant_cingulate")$anatomy, "Caudal Anterior Cingulate")
    expect_equal(decode_antspymm_label("dti_md_post_cingulate")$anatomy, "Posterior Cingulate")
    expect_equal(decode_antspymm_label("dti_fa_ant_int_cap")$anatomy, "Anterior Internal Capsule")
    expect_equal(decode_antspymm_label("dti_fa_sup_cor_rad")$anatomy, "Superior Corona Radiata")
    expect_equal(decode_antspymm_label("dti_fa_inf_cor_rad")$anatomy, "Inferior Corona Radiata")
    expect_equal(decode_antspymm_label("dti_fa_sup_cereb_ped")$anatomy, "Superior Cerebellar Peduncle")
    
    # Test asymmetry as a measurement
    res_asym <- decode_antspymm_label("nbm_asym_ant")
    expect_equal(res_asym$measurement, "Asymmetry")
    expect_equal(res_asym$anatomy, "Nucleus Basalis of Meynert Anterior")
  })
  
  test_that("Longest match is prioritized", {
    res_long <- decode_antspymm_label("T1Hier_vol_left_superior_temporal")
    expect_equal(res_long$anatomy, "Superior Temporal Gyrus")
    
    res_short <- decode_antspymm_label("T1Hier_vol_left_temporal")
    expect_equal(res_short$anatomy, "Temporal Lobe")
  })

  cat("\n✓ All tests completed successfully.\n")

} else {
  warning("The 'testthat' package is not installed. Please install it to run the validation suite.")
}
