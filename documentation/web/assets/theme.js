(function applySavedTheme() {
  const storageKey = 'femaster-docs-theme';
  const systemTheme = matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  let theme = systemTheme;

  try {
    const savedTheme = localStorage.getItem(storageKey);
    if (savedTheme === 'light' || savedTheme === 'dark') theme = savedTheme;
  } catch (_) {
    // Some browsers restrict localStorage for file:// pages. The system theme
    // remains a usable fallback in that case.
  }

  document.documentElement.setAttribute('data-theme', theme);
})();
