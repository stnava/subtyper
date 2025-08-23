library(subtyper)

# test_that("generate_idp_interpretation_prompt outputs valid prompt with plausibility", {
  # Generate prompt with examples
  prm <- generate_idp_interpretation_prompt(
    response_length = "medium",
    tone = "skeptical",
    n_examples = 5
  )
  
  # Basic structure
  expect_type(prm, "character")
  expect_true(nchar(prm) > 500)  # should be a long prompt
  
  # Check that plausibility requirement is present
  expect_match(prm, "plausibility")
  expect_match(prm, "range from 0 to 1")
  
  # Ensure JSON keys are specified correctly
  expect_match(prm, "'consistency'")
  expect_match(prm, "'justification'")
  expect_match(prm, "'plausibility'")
  
  # Extract all example JSONs from the prompt
  example_jsons <- stringr::str_extract_all(prm, "\\{[\\s\\S]*?\\}")[[1]]
  expect_gt(length(example_jsons), 0)
  
  # Parse each JSON example and validate schema
  for (ej in example_jsons) {
    parsed <- jsonlite::fromJSON(ej)
    expect_named(parsed, c("consistency", "justification", "plausibility"))
    
    # Check consistency is valid
    expect_true(parsed$consistency %in% c("low", "medium", "high"))
    
    # Justification must be a string
    expect_type(parsed$justification, "character")
    expect_true(nchar(parsed$justification) > 20)
    
    # Plausibility score is numeric between 0 and 1
    expect_type(parsed$plausibility, "double")
    expect_gte(parsed$plausibility, 0)
    expect_lte(parsed$plausibility, 1)
  }
  
  # Test without exemplars
  prm2 <- generate_idp_interpretation_prompt(
    response_length = "short",
    tone = "neutral",
    n_examples = 0
  )
  expect_false(grepl("Exemplar Cases", prm2))
 #})

toy_df <- data.frame(
    PC_Name = c("PC-1"),
    Perf.Dom = c("delayed recall"),
    IDP.1 = c("postcentral cortex"),
    IDP.2 = c("precentral cortex"),
    IDP.3 = c("fusiform cortex"),
    IDP.4 = c("uncinate fasciculus"),
    stringsAsFactors = FALSE
  )

prm0 = generate_idp_interpretation_prompt( 'long', 'critical', n_examples = 0)
prm1 = generate_idp_interpretation_prompt( 'long', 'skeptical', n_examples = 0)
#################
library(subtyper)
toy_df0=read.csv("UKB_train_july_2025_asym2params_energy_nc_mix_ica_sparseness_0.8_constraint_orthox0.001x1_optimizer_adam_llm_interp.csv")
toy_df=toy_df0[,sort(unique(c(1:2,grep("IDP", colnames(toy_df0)))))]
prm = generate_idp_interpretation_prompt( 'long', 'critical', n_examples = 6)
result0 <- assess_idp_consistency(
    toy_df,
    Perf.Dom = "Perf.Dom",
    idp_cols = colnames(toy_df)[grepl("^IDP\\.", colnames(toy_df))],
    prompt = prm, 
    backend = "groq",    # or "openrouter" depending on your setup
    api_key_env = "GROQ_API_KEY")
  

##########################################################################################
toy_df1=read.csv("UKB_train_july_2025_asym2params_energy_acc_mix_ica_sparseness_0.8_constraint_orthox0.001x1_optimizer_adam_llm_interp.csv")
toy_df1=toy_df1[,sort(unique(c(1:2,grep("IDP", colnames(toy_df1)))))]
result1 <- assess_idp_consistency(
    toy_df1,
    Perf.Dom = "Perf.Dom",
    idp_cols = colnames(toy_df1)[grepl("^IDP\\.", colnames(toy_df1))],
    prompt = prm, 
    backend = "groq",    # or "openrouter" depending on your setup
    api_key_env = "GROQ_API_KEY")
  

result0[,c("consistency","plausibility")]
result1[,c("consistency","plausibility")]

result2[,c("consistency","plausibility")]

mean(result0$plausibility)
mean(result1$plausibility)
mean(result2$plausibility)
