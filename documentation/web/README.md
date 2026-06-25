# FEMaster Web Documentation

This directory contains a static documentation website for the FEMaster input deck
and solver workflow. It does not require a build step.

Open `index.html` directly in a browser:

```text
documentation/web/index.html
```

The site is implemented as a small single-page app:

- `index.html` defines the page shell.
- `assets/styles.css` contains the complete visual system and responsive layout.
- `assets/app.js` contains navigation, search, command reference data, and page content.

Element images are reused from the LaTeX manual:

```text
documentation/pages/figures/element_schematics/
```
