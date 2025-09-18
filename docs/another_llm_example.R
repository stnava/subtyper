library(subtyper)
# Sample data frame
sample_df <- data.frame(
  Perf.Dom = c("Cognitive Flexibility", "Memory Recall"),
  IDP1 = c("Prefrontal Cortex Volume", "Hippocampal Volume"),
  IDP2 = c("White Matter Tracts", "Amygdala Activity"),
  stringsAsFactors = FALSE
)

myprompt = generate_idp_interpretation_prompt( "long",  'critical', 5 )
# Run the assessment, explicitly passing the API key environment variable name
result <- assess_idp_consistency(
  df = sample_df,
  Perf.Dom = "Perf.Dom",
  idp_cols = c("IDP1", "IDP2"),
  # prompt = "You are a neuroscientist. Evaluate if the IDPs are consistent with the domain based on known neuroscience. Respond in JSON with keys: consistency (yes/no), justification (brief explanation), plausibility (high/medium/low).",
  prompt = myprompt,
  backend = "groq",
  api_key_env = "GROQ_API_KEY3",  # Explicitly pass the key name
  model = "llama-3.3-70b-versatile",  # Explicitly set model for clarity
  verbose = TRUE,
  max_retries=0
#  cache_dir = tempdir()  # Optional: cache to a temp directory
)

# View the result
print(result)