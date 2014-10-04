
context("addlog")

test_that("addlog works", {

  a <- 50
  b <- 60
  d <- 2
  expect_equal(addlog(a,b), log(exp(a)+exp(b)))
  expect_equal(addlog(b,a), log(exp(a)+exp(b)))
  expect_equal(addlog(a,d), log(exp(a)+exp(d)))
  expect_equal(addlog(d,a), log(exp(a)+exp(d)))
  expect_equal(addlog(b,d), log(exp(b)+exp(d)))
  expect_equal(addlog(d,b), log(exp(b)+exp(d)))
  expect_equal(addlog(a,a+300), a+300);
  expect_equal(addlog(a,a-300), a);
  expect_equal(addlog(a+300,a), a+300);
  expect_equal(addlog(a-300,a), a);

})
