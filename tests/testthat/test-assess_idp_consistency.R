library(testthat)
library(stringr)
library(jsonlite)

test_that("parse_llm_response handles valid JSON", {
  cat("Testing parse_llm_response with valid JSON input...\n")
  response_content <- '```json
  {"consistency": "high", "justification": "Valid", "plausibility": 0.9}
  ```'
  required_fields <- c("consistency", "justification", "plausibility")
  cat("Parsing JSON response...\n")
  result <- parse_llm_response(response_content, required_fields)
  cat("Checking parsed result against expected output...\n")
  expect_equal(result, list(consistency = "high", justification = "Valid", plausibility = 0.9))
  cat("Test completed: valid JSON parsing successful.\n")
})

test_that("parse_llm_response handles malformed JSON", {
  cat("Testing parse_llm_response with malformed JSON...\n")
  response_content <- '```json
  {invalid_json}
  ```'
  required_fields <- c("consistency", "justification")
  cat("Parsing malformed JSON response...\n")
  result <- parse_llm_response(response_content, required_fields)
  cat("Checking result returns NA for all fields...\n")
  expect_equal(result, list(consistency = NA, justification = NA))
  cat("Test completed: malformed JSON handled correctly.\n")
})

test_that("parse_llm_response handles missing fields", {
  cat("Testing parse_llm_response with JSON missing fields...\n")
  response_content <- '```json
  {"consistency": "low"}
  ```'
  required_fields <- c("consistency", "justification", "plausibility")
  cat("Parsing JSON with missing fields...\n")
  result <- parse_llm_response(response_content, required_fields)
  cat("Checking result includes NA for missing fields...\n")
  expect_equal(result, list(consistency = "low", justification = NA, plausibility = NA))
  cat("Test completed: missing fields handled correctly.\n")
})

test_that("cache_read retrieves existing cache", {
  cat("Testing cache_read with existing cache file...\n")
  cache_dir <- tempdir()
  key <- "test_key"
  cache_file <- file.path(cache_dir, "test_key.json")
  test_data <- list(consistency = "high", justification = "Test")
  cat("Writing test data to cache file...\n")
  jsonlite::write_json(test_data, cache_file, auto_unbox = TRUE)
  cat("Reading cache file...\n")
  result <- cache_read(cache_dir, key)
  cat("Checking retrieved data matches expected...\n")
  expect_equal(result, test_data)
  cat("Test completed: cache read successful.\n")
})

test_that("cache_read returns NULL for non-existent file", {
  cat("Testing cache_read with non-existent cache file...\n")
  cache_dir <- tempdir()
  key <- "non_existent_key"
  cat("Attempting to read non-existent cache file...\n")
  result <- cache_read(cache_dir, key)
  cat("Checking result is NULL...\n")
  expect_null(result)
  cat("Test completed: non-existent file handled correctly.\n")
})

test_that("cache_write creates and writes to cache file", {
  cat("Testing cache_write functionality...\n")
  cache_dir <- tempdir()
  key <- "test_key"
  parsed <- list(consistency = "high", justification = "Test")
  cat("Writing parsed data to cache file...\n")
  cache_write(cache_dir, key, parsed)
  cache_file <- file.path(cache_dir, paste0(key, ".json"))
  cat("Checking if cache file was created...\n")
  expect_true(file.exists(cache_file))
  cat("Reading and verifying cache file content...\n")
  read_data <- jsonlite::read_json(cache_file, simplifyVector = TRUE)
  expect_equal(read_data, parsed)
  cat("Test completed: cache write successful.\n")
})

test_that("build_api_config sets correct defaults for groq", {
  cat("Testing build_api_config with groq defaults...\n")
  Sys.setenv(GROQ_API_KEY = "test_key")
  cat("Building API config for groq backend...\n")
  config <- build_api_config(backend = "groq", skip_api_key_check = TRUE)
  cat("Checking config values...\n")
  expect_equal(config$backend, "groq")
  expect_equal(config$api_url, "https://api.groq.com/openai/v1/chat/completions")
  expect_equal(config$model, "llama3-70b-8192")
  expect_equal(config$api_key, "test_key")
  cat("Test completed: groq config defaults set correctly.\n")
})

test_that("build_api_config errors without API key", {
  cat("Testing build_api_config with missing API key...\n")
  Sys.setenv(GROQ_API_KEY = "")
  cat("Attempting to build config without API key...\n")
  expect_error(
    build_api_config(backend = "groq", skip_api_key_check = FALSE),
    "API key not found"
  )
  cat("Test completed: missing API key error handled correctly.\n")
})

test_that("build_api_config uses custom model", {
  cat("Testing build_api_config with custom model...\n")
  Sys.setenv(GROQ_API_KEY = "test_key")
  cat("Building config with custom model...\n")
  config <- build_api_config(backend = "groq", model = "custom-model", skip_api_key_check = TRUE)
  cat("Checking custom model in config...\n")
  expect_equal(config$model, "custom-model")
  cat("Test completed: custom model set correctly.\n")
})

library(httr)
library(glue)
library(digest)

if ( FALSE ) {
test_that("query_api_robust handles caching", {
  cat("Testing query_api_robust with cached response...\n")
  cache_dir <- tempdir()
  key <- digest::digest("test_prompt")
  cached_data <- list(consistency = "high", justification = "Cached")
  cat("Writing cached data to file...\n")
  jsonlite::write_json(cached_data, file.path(cache_dir, paste0(key, ".json")), auto_unbox = TRUE)
  cat("Building API config...\n")
  config <- build_api_config("groq", skip_api_key_check = TRUE)
  cat("Running query_api_robust with caching...\n")
  result <- query_api_robust(
    domain = "test_domain",
    idps = c("idp1", "idp2"),
    config = config,
    system_prompt = "System prompt",
    user_prompt_template = "Domain: {domain}\nIDPs: {idps}\n{extra}",
    user_input_prompt = "Test",
    required_fields = c("consistency", "justification"),
    collapse_fn = function(x) paste(x, collapse = ", "),
    cache_dir = cache_dir,
    verbose = FALSE
  )
  cat("Checking cached result matches expected...\n")
  expect_equal(result, cached_data)
  cat("Test completed: caching in query_api_robust successful.\n")
})
}

test_that("query_api_robust handles API response", {
  cat("Testing query_api_robust with API response (skipped unless API key provided)...\n")
  skip("Requires API mocking or real API key")
  cat("Building API config for groq...\n")
  config <- build_api_config("groq")
  cat("Running query_api_robust with real API call...\n")
  result <- query_api_robust(
    domain = "cognition",
    idps = c("brain_volume", "white_matter"),
    config = config,
    system_prompt = "Assess neuroscientific consistency.",
    user_prompt_template = "Domain: {domain}\nIDPs: {idps}\n{extra}",
    user_input_prompt = "Return JSON.",
    required_fields = c("consistency", "justification", "plausibility"),
    collapse_fn = function(x) paste(x, collapse = ", "),
    cache_dir = NULL,
    verbose = FALSE
  )
  cat("Checking API response contains expected fields...\n")
  expect_true(all(names(result) %in% c("consistency", "justification", "plausibility")))
  cat("Test completed: API response handled correctly.\n")
})

library(dplyr)
library(tibble)

test_that("assess_idp_consistency processes data frame", {
  cat("Testing assess_idp_consistency with data frame...\n")
  cache_dir <- tempdir()
  df <- tibble(
    Perf.Dom = c("cognition", "memory"),
    idp1 = c("brain_volume", "gray_matter"),
    idp2 = c("white_matter", "hippocampus")
  )
  cat("Creating test data frame and pre-caching response...\n")
  key <- digest::digest(glue::glue("Domain: cognition\nIDPs: brain_volume, white_matter\nAssess."))
  jsonlite::write_json(
    list(consistency = "high", justification = "Test", plausibility = 0.8),
    file.path(cache_dir, paste0(key, ".json")),
    auto_unbox = TRUE
  )
  cat("Building API config...\n")
  config <- build_api_config("groq", skip_api_key_check = TRUE)
  cat("Running assess_idp_consistency...\n")
  result <- assess_idp_consistency(
    df = df,
    Perf.Dom = "Perf.Dom",
    idp_cols = c("idp1", "idp2"),
    backend = "groq",
    cache_dir = cache_dir,
    verbose = FALSE,
    skip_api_key_check = TRUE
  )
  cat("Checking result contains expected columns and row count...\n")
  expect_true(all(c("consistency", "justification", "plausibility") %in% names(result)))
  expect_equal(nrow(result), nrow(df))
  cat("Test completed: assess_idp_consistency processed data frame correctly.\n")
})





# Load required libraries
library(dplyr)
library(tibble)
library(jsonlite)
library(stringr)
library(httr)
library(glue)
library(digest)
library(subtyper)

# Define the system prompt (replace with your `generate_idp_interpretation_prompt` if available)
system_prompt <- paste(
  "You are a neuroscientist evaluating the consistency between imaging-derived phenotypes (IDPs)",
  "and a performance domain. For a given domain and IDPs, return a JSON object with:",
  "- consistency: 'high', 'medium', or 'low' based on neuroscientific plausibility.",
  "- justification: A brief explanation of the consistency assessment.",
  "- plausibility: A numeric score (0 to 1) reflecting the likelihood of the relationship."
)

# Create a sample data frame
df <- tibble(
  Perf.Dom = c("cognition", "memory", "attention"),
  idp1 = c("total_brain_volume", "gray_matter_volume", "frontal_cortex_thickness"),
  idp2 = c("white_matter_volume", "hippocampus_volume", "cingulate_cortex_volume")
)

# Set up caching directory
cache_dir <- tempdir()  # Use temporary directory for this example
cat("Using cache directory:", cache_dir, "\n")

# Run the query
cat("Starting query with assess_idp_consistency...\n")
result <- assess_idp_consistency(
  df = df,
  Perf.Dom = "Perf.Dom",
  idp_cols = c("idp1", "idp2"),
  prompt = system_prompt,
  backend = "groq",
  api_key_env = "GROQ_API_KEY",
  model = "llama3-70b-8192",
  required_fields = c("consistency", "justification", "plausibility"),
  max_retries = 3,
  retry_delay_base = 2,
  jitter = TRUE,
  temperature = 0.1,
  user_input_prompt = "Assess neuroscientific consistency and return JSON.",
  collapse_fn = function(x) paste(na.omit(x), collapse = ", "),
  cache_dir = cache_dir,
  verbose = TRUE,
  skip_api_key_check = FALSE  # Set to TRUE for testing without API key
)

# Display results
cat("Query completed. Results:\n")
print(result)

# Save results to a CSV file
cat("Saving results to output.csv...\n")
write.csv(result, "output.csv", row.names = FALSE)
cat("Results saved successfully.\n")