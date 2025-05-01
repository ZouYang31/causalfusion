test_that("reorder_unit_first correctly reorders rows based on unit pattern", {
  df <- data.frame(
    unit = c("control1", "treated", "control2", "control3"),
    value = c(1, 5, 2, 9)
  )

  reordered <- reorder_unit_first(df, "treated")
  print(reordered)

  # Check that rows with "treated" come first
  expect_equal(reordered$unit[1], c("treated"))

  # Ensure the total row count is unchanged
  expect_equal(nrow(reordered), 4)

  # Ensure all original rows are present
  expect_setequal(reordered$unit, df$unit)

  # Ensure order of non-matching rows is preserved after matched ones
  expect_equal(reordered$unit[2:4], c("control1", "control2", "control3"))
})
