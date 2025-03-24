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
TODO