test_that("ndFit works as expected", {

  data(EcoZUR_meta)
  data(EcoZUR_count)
  ndRes <- ndFit(W = EcoZUR_count[, 1:50], # consider only the first 50 taxa to run quickly
                 data = EcoZUR_meta,
                 A = ~ Diarrhea,
                 X = ~ sex + age_months,
                 num_crossval_folds = 2,
                 num_crossfit_folds = 2,
                 sl.lib.pi = c("SL.mean"),
                 sl.lib.m = c("SL.mean"))

  expect_true(inherits(ndRes, "ndFit"))

  # also add tests that things fail when expected (if you give the wrong inputs), and that test different
  # argument values
})
