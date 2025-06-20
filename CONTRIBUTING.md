Thanks for considering phytoclass for your latest project!

Pull requests are welcomed, please open an issue on github to get started.

## Development Workflow
* identify integration test, vignette, or testthat test to work from.
* modify & test until passing
  * edit
  * `devtools::load_all()`
  * run your test||vignette
* `git commit` changes 

## Testing
* `./tests/integration/` includes user-contributed data and `.qmd` files that use these data.
  * Running the tests is accomplished by rendering the `.qmd` files using quarto or RStudio.
* `./tests/testthat/` is where unit tests will be placed, but none exist yet.

## Pushing a Release
```
# 1. Open your package project in RStudio (clone from GitHub if needed)

# 2. Make code changes in R/ and update documentation comments

# 3. Regenerate documentation
devtools::document()

# 4. Bump version (choose "patch", "minor", or "major")
usethis::use_version("patch")

# 5. Update NEWS.md with a summary of changes

# 6. Load functions for testing
devtools::load_all() 

# 7. (Optional) Set up test infrastructure if missing
usethis::use_testthat()

# 8. Create test files for existing functions
usethis::use_test("my_function")  # Replace with actual function name

# 9. Write tests using testthat::test_that()
# e.g., test_that("my_function works", { expect_equal(my_function(2, 2), 4) })

# 10. Run all tests
devtools::test()

# 11. Check package with CRAN-level checks
devtools::check()

# 12. (Optional) Check on remote platforms
devtools::check_rhub()
devtools::check_win_devel()

# 13. Build package for CRAN submission
devtools::build()

# 14. Submit to CRAN
devtools::submit_cran() 
```
