# FEMaster Web Documentation

This directory contains a static documentation website for the FEMaster input deck
and solver workflow. It does not require a build step.

Open `index.html` directly in a browser:

```text
documentation/web/index.html
```

The site is implemented as regular static pages with a shared navigation script:

- `index.html` is the overview page.
- `pages/*.html` and `pages/dsl/**/*.html` contain one documentation page per file.
- `assets/styles.css` contains the visual system and responsive layout.
- `assets/theme.js` restores the selected color theme before each page is painted.
- `assets/command-docs.js` contains generated parser details derived from `./bin/FEMaster --document --doc-all` plus manual reference data for the primary DSL commands.
- `assets/command-explanations.js` contains the command prose adapted from `documentation/document.pdf`.
- `assets/app.js` contains shared navigation, search, command metadata, command filters, and generated element cards.

Element images are reused from the LaTeX manual:

```text
documentation/pages/figures/element_schematics/
```
