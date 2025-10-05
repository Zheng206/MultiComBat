test_that("batch_matrix returns correct structure and dimensions", {
  set.seed(20251006)
  bat_raw <- sample(c("A","B","C"), 30, replace = TRUE)
  # add an unused level to check droplevels
  bat <- factor(bat_raw, levels = c("A","B","C","Z"))

  out <- batch_matrix(bat)

  # names and basic structure
  expect_named(out, c("batch_vector","batch_matrix","batch_index","n_batches","ref"))
  expect_s3_class(out$batch_vector, "factor")
  expect_null(out$ref)

  # dimensions
  expect_equal(nrow(out$batch_matrix), length(bat))
  expect_equal(ncol(out$batch_matrix), nlevels(droplevels(bat)))

  # column names of model.matrix(~ -1 + bat) are "batA","batB","batC"
  expect_equal(colnames(out$batch_matrix), paste0("bat", levels(droplevels(bat))))

  # batch_index matches which(bat==level)
  idx_A <- which(droplevels(bat) == "A")
  expect_identical(out$batch_index[["A"]], idx_A)

  # n_batches equals lengths of indices
  expect_equal(out$n_batches[["A"]], length(idx_A))

  # one-hot rows sum to 1
  expect_true(all(rowSums(out$batch_matrix) == 1))
})

test_that("batch_matrix handles reference batch and errors", {
  set.seed(20251006)
  bat <- factor(sample(c("site1","site2","site3"), 40, TRUE))
  out_ref <- batch_matrix(bat, ref.batch = "site2")

  # ref is logical and marks site2
  expect_type(out_ref$ref, "logical")
  expect_equal(sum(out_ref$ref), sum(bat == "site2"))

  # invalid reference errors
  expect_error(batch_matrix(bat, ref.batch = "nope"),
               "Reference batch must be in the batch levels")
})

test_that("batch_matrix coerces non-factors safely", {
  set.seed(20251006)
  bat_chr <- sample(letters[1:3], 20, TRUE)
  out <- batch_matrix(bat_chr)
  expect_s3_class(out$batch_vector, "factor")
  expect_equal(ncol(out$batch_matrix), nlevels(out$batch_vector))
})
