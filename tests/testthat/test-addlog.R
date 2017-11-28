context("addlog and subtractlog")

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


test_that("subtractlog works", {

  a <- 60
  b <- 50
  d <- 2
  e <- 5
  expect_equal(subtractlog(a,b), log(exp(a)-exp(b)))
  expect_equal(subtractlog(a,d), log(exp(a)-exp(d)))
  expect_equal(subtractlog(b,d), log(exp(b)-exp(d)))
  expect_equal(subtractlog(e,d), log(exp(e)-exp(d)))
  expect_equal(subtractlog(e,d), log(exp(e)-exp(d)))
  expect_equal(subtractlog(a+300,a), 360);
  expect_equal(subtractlog(a,a-300), a);

  expect_equal(subtractlog(0.0, -a), log(1.0-exp(-a)));
  expect_equal(subtractlog(0.0, -b), log(1.0-exp(-b)));
  expect_equal(subtractlog(0.0, -d), log(1.0-exp(-d)));

})
