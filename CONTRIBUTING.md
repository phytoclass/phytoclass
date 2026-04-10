# NA

Thanks for considering phytoclass for your latest project!

Pull requests are welcomed, please open an issue on github to get
started.

## Development Workflow

- identify integration test, vignette, or testthat test to work from.
- modify & test until passing
  - edit
  - `devtools::load_all()`
  - run your test\|\|vignette
- `git commit` changes

## Testing

- `./tests/integration/` includes user-contributed data and `.qmd` files
  that use these data.
  - Running the tests is accomplished by rendering the `.qmd` files
    using quarto or RStudio.
- `./tests/testthat/` is where unit tests will be placed, but none exist
  yet.

## Updating the default F matrix

After updating the `.csv` in the `data-raw/` directory, run the
following to update the `.rda` object:

``` r
Fm <- read.csv("data-raw/Fm.csv", row.names=1)
save(Fm, file="data/Fm.rda")
```

## Pushing a Release

    # 1. Open your package project in RStudio (clone from GitHub if needed)

    # 2. Update NEWS.md with a summary of changes

    # 3. Regenerate documentation
    devtools::document()

    # 4. Check package with CRAN-level checks
    devtools::check()

    # 5. (Optional) Check on remote platforms
    devtools::check_rhub()
    devtools::check_win_devel()

    # 6. Bump version if needed
    #    Check that NEWS.md and DESCRIPTION file match.
    #    (choose "patch", "minor", or "major")
    usethis::use_version("patch")

    # 7. Submit to CRAN
    devtools::submit_cran()
