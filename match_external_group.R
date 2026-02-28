library("dplyr")
library(vecmatch)
library("MatchIt")

# Create database
generate_db <- function(metadata,
                 group_var = NULL, 
                 project_var = NULL, 
                 sex_var = NULL, 
                 age_interval = NULL) {
  
  # Remove NAs in sex and age
  metadata_vecmatch <- metadata %>%
    filter(!is.na(sex), !is.na(age))

  # Conditional filters
  if (!is.null(group_var)) {
    metadata_vecmatch <- metadata_vecmatch %>%
      filter(group_var == group_test)
  }

  if (!is.null(sex_var)) {
    metadata_vecmatch <- metadata_vecmatch %>%
      filter(sex_var == sex_var)
  }
  
  if (!is.null(project_var)) {
    metadata_vecmatch <- metadata_vecmatch %>%
      filter(project_var == project)
  }
  
  if (!is.null(age_interval) && length(age_interval) == 2) {
    metadata_vecmatch <- metadata_vecmatch %>%
      filter(age >= age_interval[1], age <= age_interval[2])
  }
  # Add the batch column (category)
  metadata_vecmatch["category"] <- "db"
  metadata_vecmatch <- metadata_vecmatch[c("sample_name","group_test", "sex", "age", "origin", "category")]
  return(metadata_vecmatch)
}


# Wrapper for running vecmatch using the metadata_db and giving the table with the expected distribution
run_vecmatch <- function(combined,
                         significance_raincloud = "t_test",
                         formula = NULL, # formula object
                         caliper = 1
) {
  formula = category ~ age * sex
  
  # 1. Calculate the Generalized propensity scores
  gps_matrix <- estimate_gps(formula = formula,
                             data = combined,
                             method = "vglm",
                             reference = reference
  )

  # 2. Define the Common Support Region (CSR)
  csr_matrix <- csregion(gps_matrix)

  # 3. Matching on the generalized propensity scores
  matched_data <- match_gps(
    csmatrix = csr_matrix,
    reference = reference,
    caliper = caliper
  )
  return(list(matched_data, csr_matrix, gps_matrix))
}


# Load metadata slim table
metadata <- read.table("metadata_slim.csv", sep=";", header=TRUE)

# Generate reference metadata table to be compared to a group defined by user
# In this case we want to select Controls from any project
# But we could select any disease present in the database, any specific project,
# specific age intervals and also specific project
metadata_db <- generate_db(metadata, 
                           group_var = "Controls") # Select based on disease group


# Test first against simulated data on patients.
# Here we define a population of 120 controls, age normally distributed with mean=60 and std=15, random
# sex sampling. 
set.seed(123)  # for reproducibility

# Parameters
n_samples <- 60          # total number of samples
age_mean <- 65             # average age
age_sd <- 10               # standard deviation for age

# Create simulated_data
simulated_data <- data.frame(
  sample_name = 1:n_samples,
  group_test = "Controls",
  sex = sample(c("F", "M"), n_samples, replace = TRUE),  # random F/M
  age = round(rnorm(n_samples, mean = age_mean, sd = age_sd)),  # normally distributed
  # age = round(runif(n_samples, min = 25, max = 90)), # randomly distributed in a given interval
  origin = "Austria",
  category = "simulated"
)

metadata_db

# Concatenate the 2 tables
merged_data <- rbind(metadata_db, simulated_data)

# Use vecmatch visualization function to show the distribution of the 2 groups
raincloud(
  data = merged_data,
  y = age,
  group = category,
  significance = "t_test",
  sig_label_color = TRUE
)


# Check if it is possible to perform 1-to-1 Sex-Age matching. Compare number of samples matched to total number of samples
m.out <- matchit(as.factor(category) ~ sex, data = merged_data,exact = ~sex+age) # we set we wanna keep the exact 1-to-1 sex and age
matched_exact_data <- match.data(m.out)
delta <- sum(merged_data$category == "simulated") - sum(matched_exact_data$category == "simulated")
# Use raincloud visualization raincloud function
raincloud(
  data = as.data.frame(matched_exact_data),
  y = age,
  group = category,
  significance = "t_test",
  sig_label_color = TRUE
)
print(delta)

# If delta is > 0, then result is suboptimal cause we wanna keep 100% of the samples in the dataset
# then we ignore this result and use vecmatch.
# Run the vecmatch package
reference <- "simulated"
formula = category ~ age * sex # This has to be defined outside the function (defining it inside actually creates some issues)
results <- run_vecmatch(combined =  merged_data)
matched_data <- results[[1]]
csr_matrix <- results[[2]]
gps_matrix <- results[[3]]


# Check distributions of matched dataset
raincloud(
  data = matched_data,
  y = age,
  group = category,
  significance = "t_test",
  sig_label_color = TRUE
)


# Check quality
# Show quality
balqual(matched_data,
        formula,
        statistic = "max"
)

# Check distributions of matched dataset


